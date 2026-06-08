# Theoretical support for the 2D-EGT / ECRAM manuscript

**Project:** `project_ecram` · **Created:** 2026-06-07
**Manuscript:** `project_ecram/ref/Draft.pdf` (Results & Discussion, EGT on 2D channels)

This note provides direct, first-principles theoretical support for the central
mechanistic claim of the manuscript and **fills the empty placeholder panel
Figure 5F** ("Simulation of interfacial ion concentration impact on the OCV of
the material").

All quantitative results come from a 1D Poisson–Nernst–Planck (PNP) solver
(`sim/pnp_core.py`) that is **validated against the analytic Gouy–Chapman–Stern
equilibrium** (reproduces the Boltzmann interfacial concentration to <0.1 %),
plus standard Nernst / Butler–Volmer electrochemistry. A COMSOL counterpart
(`sim/pulse.m`) reproduces the same accumulation in the group's existing
MATLAB+COMSOL pipeline.

---

## 1. What the manuscript needs theory for

The paper shows one EGT device operating in two modes:

| Mode | Physics | Volatility | Onset |
|------|---------|-----------|-------|
| **EDLT** | electrostatic double-layer gating | volatile | low V_G |
| **ECRAM / ECD** | Li⁺ intercalation + anion oxidation | non-volatile | high V_G |

The **headline finding** (Figs 2–4): a critical V_G separates the modes
(~4 V single-pulse for graphene, ~5 V for MoS₂), **but cumulative/sequential
pulsing lowers that threshold to ~2 V**, reproducibly, across devices and
materials.

The **proposed mechanism** (Fig 5): with closely spaced pulses the interfacial
ions do not fully relax between pulses, so **interfacial Li⁺ concentration
builds up ("pre-loads")**, and that higher local concentration lowers the
voltage required to intercalate. **Figure 5F is reserved (empty) for a
simulation of exactly this.**

---

## 2. The governing physics (one coherent chain)

**(a) Double-layer transport — PNP.** Cation (Li⁺, z=+1) and anion (z=−1) in
the electrolyte obey Nernst–Planck transport coupled to Poisson:

```
∂c_i/∂t = ∇·[ D_i ( ∇c_i + z_i (F/RT) c_i ∇φ ) ]
−∇·(ε₀ε_r ∇φ) = F Σ z_i c_i
```

Electrode at x=0 is **blocking** (no ion flux) with a **Stern layer** Robin
condition relating the metal potential φ_M(t) to the diffuse-layer potential;
x=L is a bulk reservoir. This is the same model family as the group's `ddl/`
work, here driven by a **pulse train** φ_M(t).

**(b) Two timescales set the accumulation.** The double layer charges fast,
τ_c ≈ λ_D·L/D, while salt redistributes over the bulk on τ_D ≈ L²/D ≫ τ_c.
Accumulation occurs when the **pulse period < τ_D** (and, for the interfacial
layer itself, < the EDL discharge time). 

> **Physical mapping to the experiment.** For a solid Li-salt electrolyte
> (D_Li ≈ 10⁻¹³ m²/s) over a µm-scale gate electrolyte, **τ_D ≈ tens of
> seconds**. Hence 1-s pulses at ~0.5 Hz are firmly in the accumulation regime,
> whereas the paper's *single-pulse* characterization deliberately waits **300 s
> between pulses** (Fig 2C) — i.e. ≫ τ_D — to guarantee a full reset. The model
> and the measurement protocol are consistent.

**(c) Interfacial concentration → OCV (Nernst).** The intercalation half-reaction
`Li⁺ + e⁻ + ⟨host⟩ ⇌ Li⟨host⟩` has an equilibrium (open-circuit) potential that
depends on the **local interfacial Li⁺ activity**:

```
ΔE_OC = (RT/F) ln( c_+(0) / c_∞ )      ≈ 59 mV per decade of enrichment
```

Pre-loading the interface raises c_+(0) and therefore raises the thermodynamic
driving force (OCV) for intercalation. **This is the content of Fig 5F.**

**(d) Faradaic threshold → Butler–Volmer.** The intercalation (cathodic) current
scales with the interfacial reactant concentration:

```
i_F = i₀ (c_+(0)/c_ref) [ exp(α_a F η/RT) − exp(−α_c F η/RT) ]
```

For a fixed onset current i*, the required interfacial overpotential drops by

```
Δ|η_th| = (RT/α_c F) ln( c_+(0)/c_∞ )
```

So a pre-loaded interface intercalates at lower applied voltage — the mechanism
behind the cumulative-pulse threshold reduction in Figs 3–4.

---

## 3. Results (validated, first-principles)

Run: `cd sim && python3 ecram_theory.py` → figures in `../fig/`.

### Fig 1 — `fig/fig1_accumulation.png`  (supports Figs 2B, 3, 5C–E)
- **Single pulse:** interfacial c_+(0) spikes to ~11× bulk then **fully relaxes
  to 1** → volatile EDL response.
- **Fast train (high rep rate):** peaks *and* between-pulse baseline **ratchet
  upward**; the interface never resets → **pre-loading**.
- **Slow train (low rep rate):** each pulse relaxes independently → no
  accumulation. (Direct demonstration of the rep-rate dependence the paper
  invokes.)
- Pre-load baseline rises 4.4 → **5.9× bulk** over 12 pulses (peaks 11 → 15×).

### Fig 2 — `fig/fig2_OCV_fig5F.png`  ★ **fills manuscript Fig 5F**
- Nernst mapping ΔE_OC = (RT/F)·ln(c_+/c_∞), the ~59 mV/decade line, with the
  PNP-derived pulse-train states marked.
- Driving-force gain vs pulse number: **up to ≈ +46 mV** of intercalation OCV
  purely from interfacial pre-loading.

### Fig 3 — `fig/fig3_threshold.png`  (supports Figs 3, 4)
- Butler–Volmer onset curves shift to lower overpotential as pre-load increases.
- Interfacial intercalation threshold reduction vs pulse number, anchored at the
  isolated single pulse (0 mV): reaches **≈ 91 mV at the interface** for the
  fully pre-loaded state (vs a fully-relaxed single pulse).

> **Interface vs device voltage.** The ~46–91 mV figures are *interfacial*
> (diffuse-layer / overpotential) quantities. At the device level the gate
> voltage V_G is partitioned across the Stern layer, the series electrolyte, and
> the channel; the same fractional pre-loading therefore appears as a larger V_G
> reduction, monotonic in pulse number — consistent with the volt-scale
> single→cumulative shifts in Figs 3–4. A quantitative V_G mapping needs the
> series partition (next step, §4).

---

## 4. Files & next steps

**Delivered**
- `sim/pnp_core.py` — validated 1D PNP solver (equilibrium self-test: PASS).
- `sim/ecram_theory.py` — produces Figs 1–3 above.
- `sim/explore_pulse.py` — diagnostic (rep-rate dependence).
- `sim/pulse.m` — COMSOL (LiveLink) pulse-train PNP for the group's pipeline;
  reproduces Fig 1 accumulation. Run with `N_PULSE=1` for the single-pulse
  control.

**Natural extensions (if wanted)**
1. **Coupled BV in COMSOL** — add the Butler–Volmer intercalation flux at x=0
   using the `TertiaryCurrentDistributionNernstPlanck` + `ElectrodeReaction`
   pattern from `ddl_bv/ddlBV0.m`, to get the threshold V_G(N) curve from a
   single first-principles solve.
2. **Device-level V_G mapping** — add the Stern + series + channel partition to
   convert interfacial overpotential reductions into V_G threshold reductions
   (target: reproduce 4–5 V → ~2 V).
3. **Pre-formed EDL bias** (manuscript Fig 5B,D) — hold a small DC bias, then
   pulse; the model directly predicts the enhanced response.
