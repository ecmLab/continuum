"""
Preliminary process-level digital twin for the WPI alkaline-leaching ->
HCl-precipitation Si/PV recycling sequence.

Two modules:
  (1) Si dissolution kinetics in NaOH -- surface-reaction-controlled
      shrinking-core, Arrhenius temperature dependence, first-order in [NaOH].
      Calibrated against the WPI endpoint (Wang & Ma, Resour. Conserv. Recycl.
      2024, 211, 107887): ~98 % Si conversion at 5 M NaOH, 80 C, 60 min.

  (2) pH-controlled selective precipitation -- aqueous-equilibrium
      speciation of Si(IV) and Al(III) using textbook log-K values
      for SiO2(am)/silicate and Al(OH)3(am)/aluminate equilibria.

Generates two publication-quality PDF figures into ../drafts/figs/.
"""

import os

import numpy as np
import matplotlib.pyplot as plt

# ---------- output location ----------
HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.normpath(os.path.join(HERE, "..", "drafts", "figs"))
os.makedirs(FIGDIR, exist_ok=True)


# =========================================================================
# Module 1: Si dissolution kinetics in NaOH (shrinking-core, surface rxn)
# =========================================================================
# Si + 2 NaOH + H2O -> Na2SiO3 + 2 H2
#
# Surface-reaction-controlled shrinking-core:
#   1 - (1 - X)^(1/3) = k * t
# with
#   k(T, [NaOH]) = k0 * [NaOH]^n * exp(-Ea / (R T))
#
# Calibration: WPI report ~98 % conversion at C = 5 M, T = 353 K, t = 60 min.
# Choose Ea = 60 kJ/mol (typical for Si dissolution in alkaline media),
# n = 1, and back out k0.

R_GAS = 8.314          # J / (mol K)
EA = 60.0e3            # J / mol
N_ORDER = 1.0          # order in [NaOH]


def _k(C_naoh, T):
    """Shrinking-core rate constant k(T, C) [1/min]."""
    return K0 * (C_naoh ** N_ORDER) * np.exp(-EA / (R_GAS * T))


# Calibrate k0 from WPI endpoint
_X_target = 0.98
_C_cal = 5.0           # M
_T_cal = 353.0         # K (80 C)
_t_cal = 60.0          # min
_k_cal = (1.0 - (1.0 - _X_target) ** (1.0 / 3.0)) / _t_cal
K0 = _k_cal / ((_C_cal ** N_ORDER) * np.exp(-EA / (R_GAS * _T_cal)))


def conversion(t, C_naoh, T):
    """Si conversion X(t) under surface-reaction-controlled shrinking-core."""
    kt = _k(C_naoh, T) * t
    inner = np.clip(1.0 - kt, 0.0, 1.0)
    return 1.0 - inner ** 3


def plot_leaching():
    t = np.linspace(0.0, 120.0, 400)
    T_C = 80.0
    T = T_C + 273.15
    concs = [3.0, 5.0, 7.0]
    colors = ["#1f77b4", "#d62728", "#2ca02c"]

    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    for C, col in zip(concs, colors):
        X = conversion(t, C, T)
        ax.plot(t, 100.0 * X, color=col, lw=2.0,
                label=fr"[NaOH] = {C:.0f} M")

    # WPI calibration point
    ax.plot([60.0], [98.0], marker="o", ms=10, mfc="none",
            mec="black", mew=2.0, ls="none",
            label="WPI experiment (5 M, 80 $^\\circ$C, 60 min)")

    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Si conversion (%)")
    ax.set_xlim(0.0, 120.0)
    ax.set_ylim(0.0, 102.0)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right", framealpha=0.95, fontsize=9)
    ax.set_title("Predicted Si leaching dynamics in NaOH at 80 $^\\circ$C",
                 fontsize=11)

    fig.tight_layout()
    out = os.path.join(FIGDIR, "leaching_kinetics.pdf")
    fig.savefig(out)
    plt.close(fig)
    return out


# =========================================================================
# Module 2: pH-controlled solubility of Si(IV) and Al(III)
# =========================================================================
# Goal: identify the pH window in which Si precipitates as amorphous
# SiO2 while Al(III) remains fully dissolved (selectivity for Si recovery
# without Al co-precipitation).
#
# Si(IV) total dissolved concentration (simplified):
#   For pH <= 9, controlled by amorphous-SiO2 solubility, [Si(OH)4] ~ 2 mM.
#   For pH > 9, silicate ions (H3SiO4-, H2SiO4^2-) dominate, solubility
#   rises ~ 1 decade per pH unit (pKa1 ~ 9.84, pKa2 ~ 13.1).
#
# Al(III) total dissolved concentration (amphoteric):
#   At low pH, [Al^3+] dominates and goes very high (Ksp Al(OH)3 ~ 3e-34).
#   Between pH ~4-9, controlled by Al(OH)3 solubility, ~1e-7 to 1e-5 M.
#   At high pH, Al(OH)4^- dominates, solubility rises again.
#
# These are textbook log-K solubility envelopes for amorphous hydroxides,
# rough but sufficient for a preliminary selectivity map.

pKa1_Si = 9.84
pKa2_Si = 13.10
Si_neutral = 2.0e-3        # M, amorphous SiO2 solubility ~ 2 mM

logKsp_AlOH3 = -33.5       # log [Al^3+][OH-]^3
logK4_Al = -1.0            # K4 = [Al(OH)4-]/[OH-] ~ 0.1 with Al(OH)3(s)


def si_solubility(pH):
    """Total dissolved Si(IV) (mol/L) above amorphous SiO2(s)."""
    # neutral H4SiO4 species at fixed activity
    H4SiO4 = np.full_like(pH, Si_neutral, dtype=float)
    # H3SiO4-
    H3SiO4 = H4SiO4 * 10.0 ** (pH - pKa1_Si)
    # H2SiO4^2-
    H2SiO4 = H4SiO4 * 10.0 ** (2.0 * pH - pKa1_Si - pKa2_Si)
    return H4SiO4 + H3SiO4 + H2SiO4


def al_solubility(pH):
    """Total dissolved Al(III) (mol/L) above Al(OH)3(s)."""
    pOH = 14.0 - pH
    log_OH = -pOH
    # [Al^3+] from Al(OH)3(s) + 3 H+ -> Al^3+ + 3 H2O,
    # log[Al^3+] = logKsp - 3 log[OH-]
    log_Al3 = logKsp_AlOH3 - 3.0 * log_OH
    # [Al(OH)4-]
    log_AlOH4 = logK4_Al + log_OH
    # cap below saturation (avoid runaway at extreme pH)
    return 10.0 ** log_Al3 + 10.0 ** log_AlOH4


def plot_selectivity():
    pH = np.linspace(0.0, 14.0, 561)
    Si = si_solubility(pH)
    Al = al_solubility(pH)

    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    ax.semilogy(pH, Si, color="#1f77b4", lw=2.2,
                label="Si(IV) total dissolved")
    ax.semilogy(pH, Al, color="#d62728", lw=2.2,
                label="Al(III) total dissolved")

    # selective Si-precipitation operating window
    ax.axvspan(0.0, 2.0, color="#2ca02c", alpha=0.15,
               label="Selective SiO$_2$ recovery window (pH $<$ 2)")
    ax.axhline(Si_neutral, color="#1f77b4", ls=":", lw=1.0, alpha=0.6)

    ax.set_xlabel("pH")
    ax.set_ylabel("Equilibrium solubility (mol/L)")
    ax.set_xlim(0.0, 14.0)
    ax.set_ylim(1e-9, 1e2)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="upper center", framealpha=0.95, fontsize=9)
    ax.set_title("pH-controlled selectivity: Si precipitates, Al stays dissolved",
                 fontsize=11)

    fig.tight_layout()
    out = os.path.join(FIGDIR, "selective_precipitation.pdf")
    fig.savefig(out)
    plt.close(fig)
    return out


# =========================================================================
# main
# =========================================================================

def main():
    # Diagnostic: report the calibrated rate constants
    T = 353.15
    print("Calibration check (T = 80 C):")
    for C in [3.0, 5.0, 7.0]:
        k = _k(C, T)
        # time to 98 % conversion
        t98 = (1.0 - (1.0 - 0.98) ** (1.0 / 3.0)) / k
        print(f"  [NaOH] = {C:.1f} M  ->  k = {k:.4f} 1/min,"
              f"  t(X=0.98) = {t98:.1f} min")

    # Diagnostic: solubility envelopes at key pH values
    print("\nSolubility envelopes:")
    for pH in [0.5, 1.0, 2.0, 7.0, 11.0, 13.0]:
        s_si = si_solubility(np.array([pH]))[0]
        s_al = al_solubility(np.array([pH]))[0]
        print(f"  pH = {pH:4.1f}  ->  [Si]_tot = {s_si:.2e} M,"
              f"  [Al]_tot = {s_al:.2e} M")

    f1 = plot_leaching()
    f2 = plot_selectivity()
    print("\nWrote figures:")
    print(f"  {f1}")
    print(f"  {f2}")


if __name__ == "__main__":
    main()
