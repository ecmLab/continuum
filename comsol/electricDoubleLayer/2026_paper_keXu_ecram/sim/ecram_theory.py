"""
ecram_theory.py
---------------------------------------------------------------------------
Theoretical support for the 2D-EGT / ECRAM manuscript (project_ecram).

Three first-principles results, all driven by the validated 1D PNP solver
(pnp_core.py, checked against analytic Gouy-Chapman-Stern equilibrium):

  FIG 1  Interfacial Li+ accumulation under pulsing
         (single pulse relaxes; fast train pre-loads the interface)
         -> supports manuscript Fig 2B, Fig 3, Fig 5C-E mechanism

  FIG 2  Interfacial ion concentration -> OCV / intercalation driving force
         (Nernst mapping of the PNP enrichment)
         -> FILLS the empty placeholder panel manuscript Fig 5F

  FIG 3  Butler-Volmer intercalation threshold vs pulse number
         (pre-loading lowers the voltage needed to intercalate)
         -> supports manuscript Fig 3, Fig 4 (4-5 V single -> ~2 V cumulative)

Physical mapping (solid Li-salt electrolyte, slow Li+):
  D_Li ~ 1e-13 m^2/s, electrolyte L ~ few um  =>  tau_D = L^2/D ~ tens of s.
  1-s pulses at ~0.5 Hz  =>  pulse period << tau_D  =>  accumulation regime.
  Well-separated single pulses (300-s relaxation in the paper) => full reset.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pnp_core import PNP1D

RT_F = 8.314 * 298.15 / 96485.0 * 1000.0   # thermal voltage in mV (=25.69 mV)
ALPHA = 0.5                                  # charge-transfer coefficient

# ------------------------------------------------------------------ helpers
def pulse_train(amp, ton, toff, npulse, t0=0.02, smooth=None):
    period = ton + toff
    if smooth is None:
        smooth = 0.15 * ton

    def f(t):
        if t < t0:
            return 0.0
        k = int((t - t0) // period)
        if k >= npulse:
            return 0.0
        tp = (t - t0) % period
        up = 0.5 * (1 + np.tanh(tp / smooth)) * 0.5 * (1 + np.tanh((ton - tp) / smooth))
        return amp * up
    return f


def event_times(ton, toff, npulse, t0=0.02):
    """exact end-of-ON (peak) and end-of-OFF (baseline) sampling times."""
    period = ton + toff
    peaks, troughs = [], []
    for k in range(npulse):
        peaks.append(t0 + k * period + ton)
        troughs.append(t0 + k * period + ton + 0.98 * toff)
    return np.array(peaks), np.array(troughs)


# ------------------------------------------------------------------ model
m = PNP1D(N=180, eps=0.02, delta=0.1, Dratio=2.0, beta=5.0)
AMP = -4.0          # pulse amplitude (thermal-voltage units at the metal)
TON = 0.03          # on-time   (units of L^2/D)
TOFF = 0.012        # off-time  (fast repetition, < EDL discharge time)
NPULSE = 12

print("PNP parameters:")
print(f"  eps=lambda_D/L = {m.eps},  tau_c(EDL)~eps={m.eps},  tau_D=1")
print(f"  pulse amp={AMP} (thermal V), ton={TON}, toff={TOFF}, N={NPULSE}")

# ---- single pulse (then long rest) ----
f_single = pulse_train(AMP, TON, toff=1.4, npulse=1)
t_s = np.linspace(0, 1.0, 1000)
sol_s = m.run(f_single, t_s)
cpO_s = sol_s.y[0, :]
V_s = np.array([f_single(t) for t in t_s])

# ---- fast train ----
f_fast = pulse_train(AMP, TON, TOFF, NPULSE)
peaks_t, troughs_t = event_times(TON, TOFF, NPULSE)
t_dense = np.linspace(0, (TON + TOFF) * NPULSE + 0.5 + 0.02, 1600)
t_f = np.unique(np.concatenate([t_dense, peaks_t, troughs_t]))
sol_f = m.run(f_fast, t_f)
cpO_f = sol_f.y[0, :]
V_f = np.array([f_fast(t) for t in t_f])

# enrichment per pulse number: baseline (trough) and peak interfacial conc.
def at(times):
    idx = [np.argmin(np.abs(sol_f.t - tt)) for tt in times]
    return cpO_f[idx]

R_base = at(troughs_t)     # between-pulse pre-load level vs pulse number
R_peak = at(peaks_t)
Npulse_axis = np.arange(1, NPULSE + 1)

print("\nInterfacial cation enrichment c_s/c_bulk vs pulse number (fast train):")
for n, rb, rp in zip(Npulse_axis, R_base, R_peak):
    print(f"  N={n:2d}:  baseline={rb:5.2f}   peak={rp:5.2f}")

# ---- slow train (control: full relaxation) ----
f_slow = pulse_train(AMP, TON, toff=0.12, npulse=8)
t_sl = np.linspace(0, (TON + 0.12) * 8 + 0.3, 1400)
sol_sl = m.run(f_slow, t_sl)
cpO_sl = sol_sl.y[0, :]

# ============================================================ FIGURE 1
fig1, ax = plt.subplots(1, 3, figsize=(15, 4.2))

ax[0].plot(t_s, cpO_s, 'b', lw=1.8)
a0 = ax[0].twinx(); a0.plot(t_s, V_s, 'r--', lw=1, alpha=.5); a0.set_ylabel(r"$\phi_M$ (kT/e)", color='r')
ax[0].set_title("(A) Single pulse  →  volatile (EDL)")
ax[0].set_xlabel(r"time  ($L^2/D$)"); ax[0].set_ylabel(r"interfacial $c_{+}/c_\infty$")
ax[0].axhline(1, color='k', ls=':', lw=.8)

ax[1].plot(t_f, cpO_f, 'b', lw=1.2, label="fast train (high rep rate)")
ax[1].plot(t_sl, cpO_sl, color='orange', lw=1.2, alpha=.85, label="slow train (low rep rate)")
ax[1].axhline(1, color='k', ls=':', lw=.8)
ax[1].set_title("(B) Pulse train: fast rep rate pre-loads the interface")
ax[1].set_xlabel(r"time  ($L^2/D$)"); ax[1].set_ylabel(r"interfacial $c_{+}/c_\infty$")
ax[1].legend(loc="upper right", fontsize=9)

ax[2].plot(Npulse_axis, R_base, 'o-', color='navy', label="between-pulse baseline")
ax[2].plot(Npulse_axis, R_peak, 's--', color='steelblue', label="pulse peak")
ax[2].axhline(1, color='k', ls=':', lw=.8, label="single-pulse reset level")
ax[2].set_title("(C) Interfacial pre-load vs pulse number")
ax[2].set_xlabel("pulse number  N"); ax[2].set_ylabel(r"$c_{+}(0)/c_\infty$")
ax[2].legend(fontsize=9)
fig1.suptitle("Fig 1 — Cumulative pulsing accumulates interfacial Li$^+$ (1D PNP, first-principles)", fontsize=12)
fig1.tight_layout(rect=[0, 0, 1, 0.96])
fig1.savefig("../fig/fig1_accumulation.png", dpi=130)
print("\nsaved fig1_accumulation.png")

# ============================================================ FIGURE 2  (= Fig 5F)
cs_axis = np.logspace(0, 1.2, 200)          # 1 .. ~16 x bulk
dOCV = RT_F * np.log(cs_axis)               # Nernst shift of intercalation OCV (mV)

fig2, axb = plt.subplots(1, 2, figsize=(11, 4.4))
axb[0].semilogx(cs_axis, dOCV, 'k', lw=2)
axb[0].semilogx(R_base, RT_F * np.log(R_base), 'o', color='crimson',
                label="PNP pulse-train states")
axb[0].set_xlabel(r"interfacial concentration  $c_{+}(0)/c_\infty$")
axb[0].set_ylabel(r"intercalation driving force  $\Delta E_{OC}$ (mV)")
axb[0].set_title(r"(A) Nernst: $\Delta E_{OC}=\frac{RT}{F}\ln(c_{+}/c_\infty)$  ($\approx$59 mV/decade)")
axb[0].grid(True, which='both', alpha=.3); axb[0].legend(fontsize=9)

axb[1].plot(Npulse_axis, RT_F * np.log(R_base), 'o-', color='crimson')
axb[1].set_xlabel("pulse number  N")
axb[1].set_ylabel(r"$\Delta E_{OC}$ gained (mV)")
axb[1].set_title("(B) Driving force vs pulse number (pre-loading → higher OCV)")
axb[1].grid(True, alpha=.3)
fig2.suptitle("Fig 2 — Interfacial ion concentration impact on OCV  [fills manuscript Fig 5F]", fontsize=12)
fig2.tight_layout(rect=[0, 0, 1, 0.95])
fig2.savefig("../fig/fig2_OCV.png", dpi=130)
print("saved fig2_OCV.png  (the interfacial-concentration -> OCV panel; manuscript Fig 5F)")

# ============================================================ FIGURE 3  (BV threshold)
# Butler-Volmer cathodic (intercalation) branch: i = i0 * (c_s/c_ref) * exp(alpha*f*|eta|)
# Onset (fixed target current i*) ->  |eta_th| = const - (RT/(alpha F)) ln(c_s/c_ref)
eta = np.linspace(0, 350, 400)              # interfacial overpotential magnitude (mV)
sel = [0, 3, 7, 11]                         # show N = 1,4,8,12
fig3, axc = plt.subplots(1, 2, figsize=(11, 4.4))
i_target = 1.0
for k in sel:
    cs = R_base[k]
    i = cs * np.exp(ALPHA * eta / RT_F)
    axc[0].semilogy(eta, i, lw=1.8, label=f"N={k+1}  (pre-load {cs:.1f}×)")
axc[0].axhline(i_target * (R_peak.max()), color='grey', ls='--', lw=1)
axc[0].set_xlabel(r"interfacial overpotential $|\eta|$ (mV)")
axc[0].set_ylabel("Faradaic intercalation current (a.u.)")
axc[0].set_title("(A) Butler–Volmer onset shifts left with pre-loading")
axc[0].legend(fontsize=8); axc[0].grid(True, which='both', alpha=.3)

# reduction vs the fully-relaxed single pulse (interface resets to bulk, c_s=1)
# anchor N=0 = isolated single pulse (full relaxation) -> zero reduction
N_axis2 = np.concatenate([[0], Npulse_axis])
dV_th = (RT_F / ALPHA) * np.log(np.concatenate([[1.0], R_base]))
axc[1].plot(N_axis2, dV_th, 'o-', color='darkgreen')
axc[1].annotate("isolated\nsingle pulse", (0, 0), textcoords="offset points",
                xytext=(28, 6), fontsize=8, color='grey')
axc[1].set_xlabel("pulse number  N  (0 = isolated single pulse)")
axc[1].set_ylabel(r"interfacial threshold reduction $\Delta|\eta_{th}|$ (mV)")
axc[1].set_title("(B) Intercalation threshold drops with pulse number")
axc[1].grid(True, alpha=.3)
fig3.suptitle("Fig 3 — Pre-loading lowers the Butler–Volmer intercalation threshold "
              "(mechanism behind Fig 3–4: single-pulse → cumulative)", fontsize=11)
fig3.tight_layout(rect=[0, 0, 1, 0.95])
fig3.savefig("../fig/fig3_threshold.png", dpi=130)
print("saved fig3_threshold.png")

# ------------------------------------------------------------------ summary
print("\n==== SUMMARY (numbers to quote in the theory note) ====")
print(f"  single-pulse interfacial peak           : {cpO_s.max():.1f}x bulk, relaxes to ~1")
print(f"  fast-train baseline after {NPULSE} pulses     : {R_base[-1]:.1f}x bulk (pre-load)")
print(f"  enrichment ratio baseline_N12 / bulk    : {R_base[-1]:.2f}")
print(f"  Nernst OCV gain at {R_base[-1]:.1f}x bulk         : {RT_F*np.log(R_base[-1]):.1f} mV")
print(f"  BV interfacial threshold reduction      : {(RT_F/ALPHA)*np.log(R_base[-1]/R_base[0]):.1f} mV (vs first pulse)")
print(f"  BV threshold reduction vs full reset(=1): {(RT_F/ALPHA)*np.log(R_base[-1]):.1f} mV")
