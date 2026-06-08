"""
explore_pulse.py  -- diagnostic, not a final figure.
Ask the first-principles PNP solver what happens at the electrode under
(a) a single pulse and (b) a pulse train at different repetition rates.
We watch the cation concentration at the reaction plane (OHP, node 0) and the
near-surface cation surface excess Gamma = int_0^Xc (c_p - 1) dX.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pnp_core import PNP1D


def pulse_train(amp, ton, toff, npulse, t0=0.02, smooth=None):
    """Return phiM(t): square pulses of height `amp`, baseline 0."""
    period = ton + toff
    if smooth is None:
        smooth = 0.15 * ton

    def f(t):
        if t < t0:
            return 0.0
        tp = (t - t0) % period
        k = int((t - t0) // period)
        if k >= npulse:
            return 0.0
        # smoothed square: on for [0,ton]
        up = 0.5 * (1 + np.tanh((tp) / smooth)) * 0.5 * (1 + np.tanh((ton - tp) / smooth))
        return amp * up
    return f


def surface_excess(m, cp, Xcut=0.1):
    reg = m.X <= Xcut
    return np.trapz(cp[reg] - 1.0, m.X[reg])


def cp_OHP_series(m, sol):
    return sol.y[0, :]


def run_case(m, phiM_func, t_end, npts=1200):
    t_eval = np.linspace(0, t_end, npts)
    sol = m.run(phiM_func, t_eval, rtol=1e-6, atol=1e-9)
    cpO = sol.y[0, :]
    Gam = np.array([surface_excess(m, sol.y[:m.N + 1, k]) for k in range(sol.y.shape[1])])
    phiM = np.array([phiM_func(t) for t in sol.t])
    return sol.t, phiM, cpO, Gam


if __name__ == "__main__":
    m = PNP1D(N=180, eps=0.02, delta=0.1, Dratio=2.0, beta=5.0)
    print(f"tau_c (EDL charging, dimensionless) ~ eps = {m.eps}")
    print("tau_D (bulk diffusion)              = 1.0")

    amp = -4.0
    ton = 0.03

    # single pulse, then long rest
    f_single = pulse_train(amp, ton, toff=1.5, npulse=1)
    t1, V1, cp1, G1 = run_case(m, f_single, t_end=1.2)

    # train, fast repetition (toff < tau_c): ions can't dissipate
    f_fast = pulse_train(amp, ton, toff=0.012, npulse=12)
    t2, V2, cp2, G2 = run_case(m, f_fast, t_end=1.0)

    # train, slow repetition (toff >> tau_c): full relaxation each cycle
    f_slow = pulse_train(amp, ton, toff=0.10, npulse=8)
    t3, V3, cp3, G3 = run_case(m, f_slow, t_end=1.2)

    fig, ax = plt.subplots(3, 1, figsize=(9, 9), sharex=False)
    ax[0].plot(t1, cp1, label="c_p at OHP")
    ax0b = ax[0].twinx(); ax0b.plot(t1, V1, 'r--', alpha=0.5); ax0b.set_ylabel("phiM", color='r')
    ax[0].set_title("Single pulse: interfacial cation conc. rises then RELAXES (volatile)")
    ax[0].set_ylabel("c_p(OHP)/c_bulk"); ax[0].legend(loc='upper right')

    ax[1].plot(t2, cp2, label="fast train (toff<tau_c)")
    ax1b = ax[1].twinx(); ax1b.plot(t2, V2, 'r--', alpha=0.4); ax1b.set_ylabel("phiM", color='r')
    ax[1].set_title("Fast pulse train: does the OHP concentration ratchet up?")
    ax[1].set_ylabel("c_p(OHP)/c_bulk"); ax[1].legend(loc='upper right')

    ax[2].plot(t2, G2, label="fast train")
    ax[2].plot(t3, G3, label="slow train")
    ax[2].plot(t1, G1, label="single pulse")
    ax[2].set_title("Near-surface cation surface excess Gamma(t)")
    ax[2].set_xlabel("t  (units of L^2/D)"); ax[2].set_ylabel("Gamma"); ax[2].legend()
    plt.tight_layout()
    plt.savefig("../fig/diag_pulse.png", dpi=110)
    print("saved ../fig/diag_pulse.png")

    # numeric summary: pulse-off baseline of c_p(OHP) per cycle for fast train
    print("\nfast train: c_p(OHP) min within each off-window (baseline ratchet?)")
    period = ton + 0.012
    for k in range(12):
        t0 = 0.02 + k * period + ton
        t1w = t0 + 0.012
        msk = (t2 >= t0) & (t2 <= t1w)
        if msk.any():
            print(f"  pulse {k+1:2d}: off-baseline c_p(OHP) = {cp2[msk].min():.4f},  peak = {cp2[(t2>=0.02+k*period)&(t2<=t0)].max() if ((t2>=0.02+k*period)&(t2<=t0)).any() else float('nan'):.3f}")
