# TECM Development Project Plan (One-Month Workplan)

Last updated: 2025-09-25

## Overview
- Purpose: Deliver a robust, documented TECM framework grounded in the paper’s baseline (Phase_0), and extend it in focused phases: binary species, Poisson/EDL, non‑dimensionalization, and further capabilities.
- Approach: Progress in clearly-scoped phases with specs, code, tests, and documentation; maintain exact notation/signs per the paper for baseline, and document deviations as implementation notes.
- Status summary:
  - Phase_0 (baseline understanding/doc) — improved, but still open for refinement; we may revisit.
  - Phase_1/2/3 (binary/EDL/nondim) — partial code present, not yet development‑ready.

## Timeline (approximate)
- Week 1: Phase_0 closeout; Phase_1 spec and core implementation; start Phase_2 spec.
- Week 2: Phase_2 implementation (Poisson/EDL) + validation; Phase_1 validations.
- Week 3: Phase_3 nondimensionalization integration + validations; solver guidance draft.
- Week 4: Tidy‑up; additional phases (mechanics/porous/interface breadth) scoped and queued; doc polishing and sign‑off.

---

## Phase_0 — Baseline Understanding (Paper Alignment)
- Goal: Reach complete, shared understanding of the baseline equations and conjugates and their exact mapping to code.
- Deliverables:
  - Chapter 1 (baseline) aligned with paper; conjugates (P, i, j, h) with correct signs, EN charge, and implementation notes.
  - Knowledge snapshot (see Appendix A) stored with notation, conjugates, and sign conventions.
- Tasks:
  - Finish equation review and ensure doc matches code choices and vice versa.
  - Extract a one‑page “solver posture” for baseline (scaling order, BCs, variable list).
- Acceptance:
  - You confirm Chapter 1 is correct and sufficient; no open sign/notation issues remain.

## Phase_1 — Binary Species Transport (priority)
- Scope: Extend transport from single effective species to binary (c+, c−) with (μ+, μ−), minimal complexity (no full Stefan–Maxwell).
- Equations (development form):
  - j_i = −M_i ∇μ_i + (t_i)/(z_i F) i; i = −σ^0 ∇Φ − (σ^0/F) Σ_{i∈{+,−}} t_i ∇μ_i.
- Design choices:
  - Keep it binary-only to serve EDL and a subset of liquid/soft electrolytes.
  - Preserve paper’s baseline style and conjugate definitions; document any implementation notes.
- Work items:
  - Spec: finalize variable names, units, t+ + t− = 1 constraints, and parameter inputs.
  - Materials: ensure EntropicChemicalEnergyDensity and ChemicalPotential support per‑species μ_i; verify SpeciesFluxNoEN interfaces; unify with EN baseline where Poisson is not used.
  - Interfaces: add stoichiometric coupling for reactions if needed; ensure flux continuity.
  - Tests/validation: simple 1D binary diffusion with migration off (sanity), migration on with constant field (analytical check).
- Acceptance:
  - Unit tests for j_i and i pass; parameter checks (valences, t_i sum) enforced; two benchmark inputs run to completion with expected trends.

## Phase_2 — Space‑Charge (Poisson) & EDL (priority)
- Scope: Add Poisson potential −∇·(ε ∇Ψ) = ρ_e and compatible binary closures to resolve EDL.
- Design choices:
  - Mesh guidance to resolve Debye length (h ≤ λ_D/3 near interfaces);
  - Permittivity ε(c,T) hook; charge density ρ_e computed from species.
- Work items:
  - Spec: boundary conditions and interface conditions (potential alignment, surface charge options);
  - Kernels/Materials: finalize PoissonEquation; ensure charge density and permittivity available; confirm CurrentDensityNoEN and SpeciesFluxNoEN work in coupled mode.
  - Validation benchmarks:
    - 1D planar EDL (Gouy–Chapman) profile comparison;
    - Electrochemical RAM transient charging (qualitative features: capacitance/retention trends).
- Acceptance:
  - EDL profile within tolerance to reference; stable convergence for a range of Debye lengths; ECRAM minimal case runs.

## Phase_3 — Non‑Dimensionalization (priority)
- Scope: Provide consistent characteristic scales (L0, t0, φ0, μ0, σ0, κ) and dimensionless kernels/materials.
- Work items:
  - Spec: map every variable and coefficient to dimensionless form; toggles for dimensional vs non‑dimensional runs;
  - Implement: confirm NonDimensionalParameters, NonDimensionalDiffusivity, and NonDimensionalDiffusion used consistently; add scaling for any missing terms;
  - Solver posture: document scaling order to improve convergence.
- Acceptance:
  - Dimensional and non‑dimensional runs produce consistent physics; observed solver robustness improvement (fewer nonlinear iterations) across selected cases.

## Phase_4 — Inelastic Mechanics and Damage
- Scope: Add viscoelastic/plastic models and damage/cohesive or phase‑field fracture (minimal viable set).
- Work items: small‑strain viscoelastic (Maxwell/KV) for electrodes; cohesive interface (traction‑separation) for debonding; stress‑assisted transport coupling check.
- Acceptance: simple load/unload hysteresis and interface debonding demos with transport coupling.

## Phase_5 — Interface Physics Breadth
- Scope: Generalize BV: asymmetric α_a, α_c; i0(c,T,σ); U(c,T,σ); optional Marcus kinetics; interphase/film growth and roughness evolution.
- Acceptance: at least one case shows overpotential asymmetry and film‑growth impact on performance.

## Phase_6 — Porous/Multiphase Structure
- Scope: Homogenized porosity/tortuosity; optional Darcy coupling; structure–transport feedback.
- Acceptance: composite electrode demo with porosity impact on effective transport and reactions.

## Phase_7 — Non‑Ideal Thermodynamics
- Scope: Activities (regular solution), gradient‑energy (Cahn–Hilliard) for phase separation (optional where relevant).
- Acceptance: phase‑separating electrode benchmark reproduces expected coexistence and interface width trends.

## Phase_8 — Numerical Strategy Guidance
- Scope: Provide a solver guide (scaling, line search, field‑split preconditioning, adaptive time/mesh strategies), with templates.
- Acceptance: guide enables out‑of‑the‑box stable runs for 3D examples.

## Phase_9 — Field Property Dependencies
- Scope: Provide σ(c,T,F), ε(c,T), κ(c,T), M(c,T) hooks; allow anisotropy and deformation effects.
- Acceptance: demos showing strong coupling (e.g., σ(F) in stretched solid electrolyte) and their impact on transport.

## Phase_10 — Thermal Boundary Physics
- Scope: Add convection/radiation BCs and anisotropic κ; latent heat where needed.
- Acceptance: thermal management demo with realistic boundary heat transfer.

## Phase_11 — Thermo–Electro Cross‑Effects
- Scope: Add Seebeck/Peltier/Thomson and Soret/Dufour couplings (minimal consistent set).
- Acceptance: verification on simple thermoelectric/thermogalvanic scenarios.

---

## Risks & Dependencies
- Equation fidelity: maintain strict alignment with paper notation/signs; deviations must be explicitly documented as implementation notes.
- Solver robustness: prioritize scaling and preconditioning guidance; introduce features behind toggles to avoid destabilizing baselines.
- Validation data: reference curves (EDL: Gouy–Chapman; binary migration sanity; ECRAM qualitative) collected early.

## Environment & Tools
- Build: MOOSE/libMesh/PETSc, clang/opt; scripts in repo; examples under `examples/`.
- Testing: `run_tests` harness; add focused inputs per phase; use gold outputs where appropriate.
- Docs: MooseDocs + LaTeX; keep `doc/codex/TECM_Development.tex` as the normative spec.

## Change Control
- Keep PRs phase‑scoped; add acceptance checks in CI (syntax, small benchmarks); tag milestones by phase.

---

## Appendix A — Essential Knowledge Snapshot (Phase_0)
- Total potential (paper): Π = ∫Ω e dV + ∫∂Ω γ dA − P, with e = ψ̇ + q − χ − ζ − μ ċ.
- Conjugates (aligned to paper and doc edits):
  - P = ∂ψ/∂F; i = −∂q/∂(∇Φ); j = −∂(q+ζ)/∂(∇μ); h = −T ∂χ/∂(∇T).
- EN charge conservation baseline: ∇·i = 0 (Poisson/EDL is an extension, not baseline).
- Implementation notes:
  - Heat: code uses log‑T potential H with ∂H/∂(∇ln T) = κ ∇T (equivalent to χ(∇T)).
  - Electrical kinetic potential q may depend on F: q = 1/2 ∇Φᵀ σ(F) ∇Φ.
- Code anchors (current repo): CurrentDensity, MassFlux, HeatFlux, FirstPiolaKirchhoffStress, DualChemicalEnergyDensity, BulkChargeTransport, FourierPotential, MechanicalEnergyDensity, MechanicalDeformationGradient, PoissonEquation (Case_1), SpeciesFluxNoEN, CurrentDensityNoEN, NonDimensionalParameters, NonDimensionalDiffusivity, NonDimensionalDiffusion.

---

## Immediate Next Steps (Week 1)
- Phase_0: quick pass to settle any remaining doc inconsistencies you spot.
- Phase_1: finalize binary spec + implement per‑species μ_i and j_i plumbing; add two sanity inputs; enforce t+ + t− = 1 and valence checks.
- Phase_2: finalize Poisson BCs and interface handling; stand up 1D EDL benchmark.
- Phase_3: review scaling coverage; decide default toggles; wire into two examples.

