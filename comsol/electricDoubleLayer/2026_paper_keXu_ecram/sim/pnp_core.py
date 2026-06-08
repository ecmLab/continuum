"""
pnp_core.py
-----------
A compact, self-validating 1D Poisson-Nernst-Planck (PNP) solver for a single
blocking electrode + Stern layer facing a bulk electrolyte reservoir.

This is the first-principles engine behind the ECRAM theoretical-support figures.
It is deliberately dimensionless (Bazant/Gouy-Chapman-Stern scaling) so the same
solve covers any material set; physical numbers are mapped back in the driver.

Scaling
    X   = x / L            (L = electrolyte thickness),  X in [0, 1]
    phi = (F/RT) * potential   (thermal-voltage units)
    c   = concentration / c_bulk
    t   = time * D / L^2       (1 time unit = bulk diffusion time tau_D = L^2/D)

Governing equations (cation z=+1, anion z=-1, equal D here)
    dc_p/dt =  d/dX[ dc_p/dX + c_p dphi/dX ]
    dc_m/dt =  d/dX[ dc_m/dX - c_m dphi/dX ]
    -eps^2 d2phi/dX2 = (1/2)(c_p - c_m),     eps = lambda_D / L

Boundary conditions
    X=0 (electrode): zero ion flux (blocking); Stern Robin for phi:
        dphi/dX|_0 = (phi_0 - phi_M(t)) / eta,   eta = eps_r*x_S/(eps_S*L)
    X=1 (bulk):      c_p=c_m=1, phi=0

Fluxes use the Scharfetter-Gummel (exponential-fit) scheme so the thin double
layer stays stable and positive.

Validation hook: at steady state under a held bias the solver must reproduce the
analytic Boltzmann/GCS profile c_p(X) = exp(-phi(X)).  test_equilibrium() checks it.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_banded


# ----------------------------------------------------------------------
# graded mesh, clustered at the electrode (X=0)
# ----------------------------------------------------------------------
def make_mesh(N, beta=4.0):
    """N+1 nodes on [0,1], refined near X=0 via a power-law stretch."""
    s = np.linspace(0.0, 1.0, N + 1)
    X = (np.expm1(beta * s) / np.expm1(beta))      # smooth clustering at 0
    X[0], X[-1] = 0.0, 1.0
    return X


def _bernoulli(t):
    """B(t) = t/(exp(t)-1), numerically safe; B(0)=1."""
    out = np.ones_like(t)
    small = np.abs(t) < 1e-10
    big = ~small
    out[big] = t[big] / np.expm1(t[big])
    return out


class PNP1D:
    def __init__(self, N=160, eps=0.02, eta=None, delta=0.1, Dratio=1.0, beta=4.0):
        """
        eps    = lambda_D / L  (double-layer thinness)
        delta  = eps_r*x_S/(eps_S*lambda_D) (Stern-to-Debye ratio); eta=delta*eps
        Dratio = D_anion / D_cation
        """
        self.X = make_mesh(N, beta)
        self.N = N
        self.eps = eps
        self.eta = delta * eps if eta is None else eta
        self.Dr = Dratio

        X = self.X
        self.d = np.diff(X)                          # node spacings, len N
        # dual-cell widths h_j around each node
        h = np.zeros(N + 1)
        h[1:-1] = 0.5 * (X[2:] - X[:-2])
        h[0] = 0.5 * (X[1] - X[0])
        h[-1] = 0.5 * (X[-1] - X[-2])
        self.h = h

    # ---- Poisson solve: phi at all nodes given (c_p - c_m) and phi_M ----
    def solve_phi(self, rho, phiM):
        """rho = c_p - c_m at nodes. Returns phi at nodes (len N+1)."""
        N, eps, eta, d, h = self.N, self.eps, self.eta, self.d, self.h
        e2 = eps * eps
        # tridiagonal system A phi = b  (banded)
        lo = np.zeros(N + 1)   # sub-diagonal
        di = np.zeros(N + 1)   # diagonal
        up = np.zeros(N + 1)   # super-diagonal
        b = np.zeros(N + 1)

        # interior nodes 1..N-1 : -eps^2[(phi_{j+1}-phi_j)/d_j - (phi_j-phi_{j-1})/d_{j-1}] = 0.5*rho*h
        for j in range(1, N):
            aR = e2 / d[j]
            aL = e2 / d[j - 1]
            lo[j] = -aL
            up[j] = -aR
            di[j] = aL + aR
            b[j] = 0.5 * rho[j] * h[j]

        # node 0 (Stern Robin):
        # -eps^2[(phi_1-phi_0)/d_0 - phi'(0)] = 0.5*rho_0*h_0 ; phi'(0)=(phi_0-phiM)/eta
        aR = e2 / d[0]
        di[0] = aR + e2 / eta
        up[0] = -aR
        b[0] = 0.5 * rho[0] * h[0] + (e2 / eta) * phiM

        # node N : Dirichlet phi=0
        di[N] = 1.0
        lo[N] = 0.0
        b[N] = 0.0

        ab = np.zeros((3, N + 1))
        ab[0, 1:] = up[:-1]      # super-diagonal
        ab[1, :] = di            # diagonal
        ab[2, :-1] = lo[1:]      # sub-diagonal
        phi = solve_banded((1, 1), ab, b)
        return phi

    # ---- right-hand side of the ODE system ----
    def rhs(self, t, y, phiM_func):
        N, d, h, Dr = self.N, self.d, self.h, self.Dr
        cp = y[:N + 1]
        cm = y[N + 1:]
        phiM = phiM_func(t)
        phi = self.solve_phi(cp - cm, phiM)

        dphi = phi[1:] - phi[:-1]                    # len N
        Bp = _bernoulli(dphi)
        Bm = _bernoulli(-dphi)

        # SG face fluxes between node j and j+1 (flux in +X direction)
        # cation (z=+1):  J = -(1/d)[ c_{j+1} B(-dphi) - c_j B(dphi) ]
        Jp = -(cp[1:] * Bm - cp[:-1] * Bp) / d
        # anion (z=-1):   J = -(1/d)[ c_{j+1} B(dphi) - c_j B(-dphi) ] * Dr
        Jm = -Dr * (cm[1:] * Bp - cm[:-1] * Bm) / d

        dcp = np.zeros(N + 1)
        dcm = np.zeros(N + 1)
        # interior nodes: h dc/dt = J_{j-1/2} - J_{j+1/2}
        dcp[1:N] = (Jp[:-1] - Jp[1:]) / h[1:N]
        dcm[1:N] = (Jm[:-1] - Jm[1:]) / h[1:N]
        # node 0: blocking (left flux=0): h0 dc/dt = -J_{1/2}
        dcp[0] = -Jp[0] / h[0]
        dcm[0] = -Jm[0] / h[0]
        # node N: Dirichlet -> frozen
        dcp[N] = 0.0
        dcm[N] = 0.0
        return np.concatenate([dcp, dcm])

    def run(self, phiM_func, t_eval, c0=None, rtol=1e-6, atol=1e-9):
        N = self.N
        if c0 is None:
            y0 = np.ones(2 * (N + 1))
        else:
            y0 = c0
        sol = solve_ivp(self.rhs, (t_eval[0], t_eval[-1]), y0,
                        t_eval=t_eval, method="BDF",
                        args=(phiM_func,), rtol=rtol, atol=atol)
        return sol


# ----------------------------------------------------------------------
# validation: steady bias must reproduce Boltzmann/GCS equilibrium
# ----------------------------------------------------------------------
def test_equilibrium():
    print("=== PNP equilibrium validation (vs analytic Boltzmann) ===")
    m = PNP1D(N=200, eps=0.02, delta=0.1, Dratio=1.0, beta=5.0)
    phiM = -3.0  # held bias, thermal-voltage units
    # integrate well past the charging time tau_c ~ eps (dimensionless)
    t_eval = np.linspace(0, 2.0, 60)
    sol = m.run(lambda t: phiM, t_eval)
    cp = sol.y[:m.N + 1, -1]
    cm = sol.y[m.N + 1:, -1]
    phi = m.solve_phi(cp - cm, phiM)

    cp_boltz = np.exp(-phi)
    cm_boltz = np.exp(phi)
    # compare on the diffuse-layer region (first ~20% of domain)
    reg = m.X < 0.2
    err_p = np.max(np.abs(cp[reg] - cp_boltz[reg]) / cp_boltz[reg])
    err_m = np.max(np.abs(cm[reg] - cm_boltz[reg]) / cm_boltz[reg])
    phi0 = phi[0]
    print(f"  held phi_M           = {phiM:.3f}  (thermal volts)")
    print(f"  diffuse-layer phi_0  = {phi0:.4f}")
    print(f"  Stern drop phi_M-phi0= {phiM - phi0:.4f}")
    print(f"  surface c_p(0)       = {cp[0]:.4f}   Boltzmann exp(-phi0)={np.exp(-phi0):.4f}")
    print(f"  max rel err cation   = {err_p:.2e}")
    print(f"  max rel err anion    = {err_m:.2e}")
    ok = err_p < 2e-2 and err_m < 2e-2
    print(f"  RESULT: {'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    test_equilibrium()
