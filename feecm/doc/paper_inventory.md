# JMPS Paper Inventory & Punch List

**Source:** `projects/ecmTest/doc/growth_deformation_kinematics.tex` (536 lines, 23 KB)
**Title (current):** *A Unified Framework for the Simulation of Solid-State Batteries*
**Target venue:** Journal of the Mechanics and Physics of Solids (JMPS)
**Scope:** Solid-state Li-metal SSB only
**Generated:** 2026-05-09

This is a Phase-1 inventory — what's drafted, what's missing, and which `projects/ecmTest/benchmarks/<sub>/` files supply each missing piece.

---

## Section-by-section status

### §1 Introduction
**Status:** 🔴 **Stub** — three bullet points, the third even cut off mid-sentence (line 44).

**What's missing:**
- Full motivation paragraph (why solid-state, what's the open problem)
- Literature review (existing simulation tools and their limitations)
- Statement of contribution (what the unified framework does that others don't)
- Paper outline ("In §2 we develop the kinematics...")
- Citations — currently zero in the document

**Action:** Full drafting required. Allocate ~2 pages JMPS format.

---

### §2 Methodology

#### §2.1 Kinematics of Deformation
**Status:** ✅ **Complete** (lines 47–99). Standard `F = F^e F^p F^g` decomposition with Mandel/Gurtin notation.

**Possible polish:** Add citation for the multiplicative decomposition lineage (Lee 1969; Rodriguez et al. 1994 for growth; Gurtin & Anand 2005). One missing label inline.

---

#### §2.2 Growth Deformation in Interphase Layer
**Status:** ✅ **Complete** (lines 101–188). Defines `F^g`, growth-weight parameters `α_i`, derives `\dot{F}^g = D^g F^g`, `S^g = Σ α_i (m_Ri ⊗ m_Ri)`.

**Possible polish:** Add a small figure schematic of the reference / intermediate / current configurations (a visual aid is conventional in JMPS kinematics sections).

---

#### §2.3 Non-fixed reference vector (incl. Axial / Areal / Swelling sub-subsections)
**Status:** ✅ **Complete** (lines 191–292). Correctly handles principal-vector transport, with a "Warning" callout about needing the incremental form when the reference vector isn't fixed.

**Possible polish:** Schematic figure showing axial vs areal growth would help reviewers visualize.

---

#### §2.4 Li Diffusion in Interphase Layer
**Status:** ✅ **Mostly complete** (lines 295–373). Anisotropic Fick's law, chemical-potential expressions, simplified flux equation.

**Issue:** Equation (line 372–373) `∂c/∂t + ∇·D∇c + D∇²c = 0` looks dimensionally inconsistent — the third term is a duplicate of the second under a divergence; check derivation.

**Action:** Verify the diffusion PDE and make sure the simplification step from lines 366→372 is correct. The stronger result is `∂c/∂t = ∇·(D∇c)`.

---

#### §2.5 Li⁺ migration in the Solid Electrolyte
**Status:** 🔴 **Heading only** (line 375). NO CONTENT.

**What's missing:**
- Constitutive equation for Li⁺ flux in SE (Nernst-Planck or simplified migration form)
- Coupling between Li⁺ in SE and Li in interphase (interfacial flux conservation)
- Charge-conservation / electroneutrality consideration
- Boundary conditions at the SE/anode interface

**Project source for figures:** `benchmarks/benchmark3_growthDiffusion/` (Li-ion migration coupling) + `benchmarks/full_cell/`

**Action:** Full drafting required. Allocate ~1 page.

---

#### §2.6 Contact problems during metal growth
**Status:** ✅ **Complete** (lines 377–465). Strong + weak forms, mortar-FEM Lagrange-multiplier formulation, conductivity-gap analogy.

**Possible polish:** Replace the bare URL on line 512 (`https://www.overleaf.com/...`) which is unrelated stray text leaked into the §3.1 caption.

---

#### §2.7 Free Energies and chemical potential
**Status:** ✅ **Complete** (lines 467–492). Helmholtz free energy decomposition, isotropic elastic-energy contribution, total chemical potential of Li.

**Possible polish:** "The derivation of the last term is done elsewhere and is not shown here" (line 492) — for JMPS this should be either expanded in an appendix or have a clear citation.

---

### §3 Benchmark Tests
**Section intro is good (line 495)** — describes the step-by-step verification approach.

#### §3.1 Plastic/creep behavior of Li metal material
**Status:** 🟡 **Figures present, commentary missing.**

**Existing figures:**
- `StrainRate.png` — true stress-strain at 298 K, four strain rates 4×10⁻⁵ to 2×10⁻² s⁻¹
- `temperature.png` — true stress-strain at strain rate 4×10⁻⁵ s⁻¹, multiple temperatures

**Stray text on line 512:** A leftover Overleaf URL embedded in the figure caption. **DELETE.**

**Project source:** `benchmarks/benchmark1_metalCreep/` (4 inputs: `compression_constTrueStrainRate.i`, `tension_constTrueStrainRate.i`, `tension_microns.i`, `R_Run.i`)

**What's missing:**
- 1–2 paragraphs interpreting the figures: what's the experimental reference, what does the agreement tell us about the constitutive model, what's the temperature range covered.
- Table of fitted Anand model parameters (typical for JMPS).
- Citation to Anand & Su (or wherever the constitutive law comes from).

**Action:** Write commentary. Possibly regenerate figures with consistent styling (axis labels, legend).

---

#### §3.2 Free growth of Li metal in the interphase layer
**Status:** 🔴 **Heading + 5 empty sub-subsections.**

| Sub-subsection | Project source | Specific input | What figure |
|---|---|---|---|
| Uniform normal growth on flat surface | `benchmarks/benchmark2_freeGrowth/` | `01_flat_uniformSwell.i` | Stress/strain vs time, deformed mesh snapshots |
| Uniform normal growth on inclined and cosine-shape surface | `benchmarks/benchmark2_freeGrowth/` | `inclined.i`, `cosine.i` | Surface displacement profiles, swelling distribution |
| Effect of initial thickness of the interlayer | `benchmarks/benchmark2_freeGrowth/` (parametric) | `swell.i` (parametric study) | Stress/displacement vs interlayer thickness |
| Li diffusion-induced growth on inclined and cosine-shape surface | `benchmarks/benchmark2_freeGrowth/` | `02_flat_diffusionSwell.i`, `diffusion.i` | Concentration field + induced deformation |
| Effect of Li diffusivity in the interlayer | `benchmarks/benchmark2_freeGrowth/` (parametric) | `02_flat_diffusionSwell_Ash.i` plus diffusivity sweep | Steady-state profile vs D |

**Note:** The "Effect of initial thickness" and "Effect of Li diffusivity" subsections are parametric studies and likely need a small wrapper script (sweep.py-style) to run inputs at multiple parameter values and post-process.

**Action:** For each subsection, run input(s), generate figure(s), write 1-paragraph commentary.

---

#### §3.3 Confined growth of Li metal between anode and SE
**Status:** 🔴 **Heading only**, no content.

**Project source:** `benchmarks/benchmark3_confinedGrowth/` (5 inputs all PASS): `01_conformalContact_mech.i`, `02_conformalContact_swell.i`, `test_circle.i`, `test_cosine.i`, `test_incline.i`

**What figure(s):**
- Mechanics-only baseline (`01_conformalContact_mech.i`): contact pressure distribution
- Confined swelling (`02_conformalContact_swell.i`): coupled mechanical+growth response
- Geometry sensitivity (the three `test_*` shapes): pressure/stress vs interface geometry

**Action:** Run, plot, write ~half page of commentary.

---

#### §3.4 Full coupling with Li-ion migration in the SE
**Status:** 🔴 **Heading only**, no content.

**Project source:**
- `benchmarks/benchmark3_growthDiffusion/benchmark3_growthDiffusion.i` (single input)
- `benchmarks/benchmark4_fullCoupled/full_model.i`, `full_model_flat.i`
- (Optionally) `benchmarks/benchmark4_contactGrowthDiffusion/` if any inputs are added

**What figure(s):**
- Steady-state Li-ion migration profile in SE
- Coupled Li metal growth + Li⁺ flux at the SE/Li interface
- Comparison with `benchmark4_contactGrowthDiffusion` if contact is also active

**Action:** Run, plot, write ~half page of commentary.

---

### §4 Li Metal Plating/Stripping under Different Conditions
**Status:** 🔴 **Empty heading on line 532, end of document.**

**Project sources (all in `benchmarks/`):**
- `plating_and_stripping/` (4 inputs): direct content
- `agLi/` (8 inputs): Ag interlayer enabling Li plating
- `dendrite_growth/` (4 inputs): extreme plating instabilities
- `mcmeeking_paper_results/` (13 inputs): literature validation against McMeeking SSB papers

**Suggested sub-structure for §4 (my proposal — your call):**
- §4.1 Galvanostatic plating with bonded contact — `plating_and_stripping/interlayer_left_bond_pore_left_bond.i`
- §4.2 Plating with porous interlayer — `plating_and_stripping/interlayer_left_bond_pore_left_bond_strip.i` (and equilibration variant)
- §4.3 Stripping under tension — `agLi/full_cell_210522/test*.i`
- §4.4 Validation against McMeeking experiments — `mcmeeking_paper_results/test_w_Ag.i`, `test_w_cathode.i`, `test_w_pore_interlayer_w_block.i`
- §4.5 Dendrite onset / growth — `dendrite_growth/DendriteGrowth.i`, `phaseField.i`

This is the **biggest writing effort** in the paper. Allocate 4–6 pages.

---

### §5 Conclusions / Discussion (currently absent)
**Status:** 🔴 **Section doesn't exist.**

**Action:** Write ~half page summarizing the framework's capabilities, limitations, and suggested next steps. JMPS expects this.

---

### Bibliography (currently absent)
**Status:** 🔴 **Zero citations in the document.**

**Action:** Build a `.bib` file. Key citations to include (rough guess from the content):
- Anand, L. — viscoplastic Li constitutive (relevant for §3.1)
- Gurtin, M.E. & Anand, L. — multiplicative decomposition with growth (for §2.1–§2.2)
- Rodriguez, E.K., Hoger, A., McCulloch, A.D. (1994) — biological growth → analogous framework
- Bucci et al., Bistri et al. — solid-state Li-metal SSB modeling
- McMeeking, R., Purkayastha — solid-state Li-metal mechanics (for §4 validation)
- Newman, J. — battery electrochemistry foundational
- Permann, C. et al. — MOOSE framework reference
- Kim, Kim, Suzuki — phase-field method (for §3 numerical-verification subsection if added)

---

## Cross-cutting actions

### Cleanup
- [ ] Remove stray Overleaf URL on line 512 (in §3.1 figure caption).
- [ ] Decide title: keep current "A Unified Framework..." or sharpen for JMPS (e.g., "A Coupled Electro-Chemo-Mechanical Framework for Solid-State Lithium-Metal Batteries").

### Numerical-methods description
JMPS reviewers may want a brief subsection on the MOOSE finite-element implementation:
- Mesh / element types used
- Time-stepping scheme  
- Solver / preconditioner
- Adaptive mesh refinement (if used)
- Mortar-FEM contact (already discussed in §2.6)

Suggested location: **new §2.8 Numerical Implementation** (1 page) OR move to an appendix.

### Software citation
For JMPS reproducibility expectations: cite the `ecm_test` app (with the consolidated post-2026-05-09 codebase). Could prepare a Zenodo DOI for the integrate-ashok branch or a release tag.

---

## Recommended Phase-2 sequencing

If we're heading toward a submitted manuscript, I'd suggest tackling sections in this order (worst-case-first to surface science blockers):

1. **§4 Plating/Stripping** — biggest unknown, may surface input-file issues that gate the rest.
2. **§3.4 Full coupling** — the multi-physics culmination figure.
3. **§3.2 Free growth** — 5 sub-subsections, mostly mechanical.
4. **§3.3 Confined growth** — short, contact-mechanics oriented.
5. **§3.1 Plastic/creep** — figures already exist; just write commentary.
6. **§2.5 Li⁺ migration in SE** — equations only, ~1 page.
7. **§1 Introduction** — write last (refer to actual content).
8. **§5 Conclusions** — write last.
9. **Bibliography** — accumulate as we go; finalize at end.

For each §3/§4 subsection the workflow is:
- Pick canonical input file(s) → run with `ecm-opt` → produce figure (matplotlib from CSV or Paraview from Exodus) → drop into doc folder → write 1-paragraph commentary.

---

## Open questions for you

1. **Title:** Keep "A Unified Framework for the Simulation of Solid-State Batteries", or sharpen mechanics emphasis (e.g., add "Electro-Chemo-Mechanical" or "Coupled")?
2. **§2.5 (Li⁺ migration in SE):** Do you have notes for this section elsewhere, or should I draft from a typical Nernst-Planck formulation?
3. **§4 sub-structure:** Agree with my proposed 4.1–4.5 split, or do you have a different organization in mind?
4. **Numerical-methods subsection:** §2.8 or appendix?
5. **Reference manager:** Do you have an existing `.bib` file with battery/mechanics citations from other work? Pointing me at it accelerates the bibliography step.
