# Paper 2: Microsphere Separation Experiments

**Title:** High-Throughput Selective Separation of Polystyrene Microspheres Using Large-Scale Aqueous Quadrupole Ion Traps

**Target Journals:**
1. Lab on a Chip (primary)
2. Analytical Chemistry
3. Electrophoresis
4. Small

## Paper Overview

This paper presents **experimental demonstration** of microsphere separation, including:

- 2 μm vs 10 μm carboxylate-modified polystyrene microspheres
- Analytical predictions + COMSOL validation + experimental results
- Device design and fabrication
- V = 2.5 V preliminary results (completed)
- V = 5-20 V experiments (IN PROGRESS - CRITICAL)
- Separation efficiency, purity, and throughput measurements

**Status:** 40% complete - REQUIRES higher voltage experiments before submission

## Files

- `paper2_microspheres.tex` - Main LaTeX manuscript
- `references.bib` - Bibliography file (shared with Paper 1)
- `build/` - Build artifacts directory (LaTeX auxiliary files)

## Compilation Instructions

### Standard LaTeX compilation:

```bash
pdflatex -output-directory=build paper2_microspheres.tex
bibtex build/paper2_microspheres
pdflatex -output-directory=build paper2_microspheres.tex
pdflatex -output-directory=build paper2_microspheres.tex
mv build/paper2_microspheres.pdf ./
```

### Using latexmk (recommended):

```bash
latexmk -pdf -outdir=build paper2_microspheres.tex
mv build/paper2_microspheres.pdf ./
```

### Cleanup:

```bash
# All build artifacts are in build/ folder
rm -rf build/*
```

## CRITICAL: Required Experiments

### ⚠️ MINIMUM VIABLE EXPERIMENTS (3 months):

1. **Higher voltage power supply** (5-20 V capability)
   - Analog Discovery 2 limited to 2.5 V
   - Need function generator or amplifier

2. **V = 5 V experiments** - MOST CRITICAL
   - 10 μm only → verify trapping
   - 2 μm only → verify rejection
   - Mixture (1:1) → demonstrate selective separation
   - **n = 5 trials minimum** for statistics

3. **Quantitative tracking**:
   - Particle counting (inlet vs trapped vs outlet)
   - Video analysis of trajectories
   - Separation efficiency: η = N_trapped / N_total
   - Collection purity: P = N_10μm / N_trapped

4. **Control experiments**:
   - V = 0 (no field) → gravity settling baseline
   - DC only (V=0, U≠0) → verify not simple electrophoresis
   - AC only (V≠0, U=0) → QIT mechanism

### 🎯 IDEAL EXPERIMENTS (6 months):

5. **Voltage sweep**: V = 2.5, 5, 10, 15, 20 V
6. **Efficiency vs voltage curve**
7. **Convergence time vs voltage**
8. **Flow rate optimization**
9. **Continuous operation** (hours-long stability)
10. **Surface charge validation** (zeta potential measurements)

### ❌ WITHOUT V ≥ 5V DATA, THIS PAPER IS WEAK

Current V = 2.5 V results are preliminary only - not sufficient for publication.

## Required Figures (To Be Created)

1. **Figure 1**: Device photos and CAD
   - (a) CAD cross-section
   - (b-c) Assembled device photos

2. **Figure 2**: Stability diagrams
   - 2 μm variants at V = 2.5, 5, 20 V
   - 10 μm variants at V = 2.5, 5, 20 V

3. **Figure 3**: COMSOL trajectories
   - Three voltages × two sizes = 6 panels
   - Showing convergence/divergence

4. **Figure 4**: Micromotion visualization
   - Zoomed trajectory showing oscillations

5. **Figure 5**: Experimental time-lapse (V = 2.5 V)
   - Current data (completed)

6. **Figure 6**: Experimental results (V = 5-20 V)
   - **NEEDS DATA** from ongoing experiments
   - Quantitative separation metrics

7. **Figure 7**: Separation efficiency vs voltage
   - **NEEDS DATA**

8. **Figure 8**: Mixture separation results
   - Size distribution before/after
   - **NEEDS DATA**

## Equipment Needed for Experiments

### Critical:
- [ ] High-voltage function generator (5-20 V, 10 kHz)
  - Options: Agilent 33220A, Tektronix AFG3000, or amplifier circuit
  - Estimated cost: $1-5K

### Helpful:
- [ ] Better imaging system for particle tracking
  - Current: Dino-Lite microscope (adequate but manual)
  - Upgrade: Automated particle counter
  - Estimated cost: $2-10K

### Supplies:
- [ ] More microspheres (multiple trials will consume samples)
  - 2 μm and 10 μm from Magsphere
  - Estimated cost: $500

**Total budget estimate: $1.5-15K depending on equipment choices**

## Key Results Expected

From analytical predictions (Table in manuscript):

| Voltage | 2 μm Status | 10 μm Status | Separation? |
|---------|-------------|--------------|-------------|
| 2.5 V | Unstable | Barely stable | ❌ No (weak) |
| 5 V | Unstable | **Stable** | ✅ **YES** |
| 10 V | Unstable | **Stable** | ✅ YES |
| 20 V | Unstable | **Deeply stable** | ✅ YES (fast) |

**SI = 2.3 predicts wide separation window (5-50 V)**

## Publication Timeline

### Option A: Sequential (RECOMMENDED)

- **Months 1-3**: Submit Paper 1 (theory) - ready now
- **Months 4-9**: Complete V = 5-20 V experiments
- **Months 10-12**: Write and submit Paper 2

### Option B: Wait and submit together

- **Months 1-6**: Theory + experiments in parallel
- **Months 6-7**: Submit both papers simultaneously

**Recommendation: Option A (submit Paper 1 first)**
- Paper 1 establishes framework (gets cited by Paper 2)
- Don't delay Paper 1 waiting for experiments
- 1 certain publication > 2 uncertain

## Key Contributions

1. ✅ First experimental demonstration of aqueous QIT microsphere separation
2. ⏳ Quantitative validation of analytical predictions (pending V ≥ 5V data)
3. ✅ Proof-of-concept device design
4. ⏳ Separation efficiency >95% at V ≥ 5V (predicted, needs validation)
5. ✅ Label-free, membrane-free alternative to FACS/DEP
6. ✅ Validates Separability Index (SI) metric from Paper 1

## Applications

- Microsphere purification (drug delivery, diagnostics)
- Future: Label-free cell sorting
- Future: Protein-conjugated microsphere separation (immunoassays)
- Future: Charge-based quality control

## Dependencies

- **Cites Paper 1** for theoretical framework
- Should reference "recently published" or "companion paper"
- Paper 1 should be published or accepted before Paper 2 submission

## Notes

- V = 2.5 V data is consistent with predictions but not conclusive
- MUST complete V ≥ 5V experiments for strong paper
- Control experiments (no field, DC only) are critical to prove mechanism
- Without quantitative data, this is just a preliminary report
