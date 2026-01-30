# Paper 1: Theoretical Framework

**Title:** Design Principles for Selective Particle Separation in Aqueous Quadrupole Ion Traps: A Unified Stability Framework

**Target Journals:**
1. Physical Review Applied (primary)
2. Separation and Purification Technology
3. Journal of Physical Chemistry Letters
4. Analytical Chemistry

## Paper Overview

This paper presents the **general theoretical framework** for aqueous QIT separation, covering:

- Damped Mathieu stability theory
- Separability Index (SI) metric
- Comprehensive analysis from ions to microspheres (5 orders of magnitude)
- Feasibility maps and design guidelines
- Negative result: Li+/Na+ separation infeasible at large scales
- Positive result: Microsphere separation highly feasible

**Status:** Ready for manuscript preparation (90% complete from thesis content)

## Files

- `paper1_framework.tex` - Main LaTeX manuscript
- `references.bib` - Bibliography file
- `build/` - Build artifacts directory (LaTeX auxiliary files)

## Compilation Instructions

### Standard LaTeX compilation:

```bash
pdflatex -output-directory=build paper1_framework.tex
bibtex build/paper1_framework
pdflatex -output-directory=build paper1_framework.tex
pdflatex -output-directory=build paper1_framework.tex
mv build/paper1_framework.pdf ./
```

### Using latexmk (recommended):

```bash
latexmk -pdf -outdir=build paper1_framework.tex
mv build/paper1_framework.pdf ./
```

### Cleanup:

```bash
# All build artifacts are in build/ folder
# To clean up:
rm -rf build/*
```

## Required Figures (To Be Created)

The manuscript includes placeholders for the following figures:

1. **Figure 1**: Conceptual schematic
   - QIT electrode geometry
   - Particle trajectories (stable vs unstable)
   - Electric potential visualization

2. **Figure 2**: Stability theory overview
   - Mathieu stability diagram (vacuum vs aqueous)
   - Effect of damping on stability boundaries

3. **Figure 3**: COMSOL validation
   - Three test cases with trajectory comparisons

4. **Figure 4**: Particle size spectrum analysis
   - Stability diagrams for ions, nanoparticles, microspheres

5. **Figure 5**: **CRITICAL** - Feasibility map
   - Particle radius vs required voltage
   - Constraint lines (breakdown, electrolysis)
   - Feasible/infeasible regions

6. **Figure 6**: Separability Index matrix
   - Heatmap showing SI for particle pairs

7. **Figure 7**: Design parameter nomograph
   - Practical parameter selection tool

8. **Figure 8**: Microsphere case study
   - 2 μm vs 10 μm at multiple voltages

## Additional Content Needed

1. **Literature update**: 2024-2025 recent papers
2. **Particle database table**: Comprehensive properties (Supplementary Material)
3. **MATLAB code**: Stability boundary calculations (Supplementary)
4. **COMSOL template**: Simulation files (Supplementary)
5. **Excel calculator**: Quick parameter estimation (Supplementary)

## Key Contributions

1. ✅ First unified stability theory for aqueous QITs
2. ✅ Separability Index (SI) - quantitative metric
3. ✅ Comprehensive feasibility analysis (ions to microspheres)
4. ✅ Scaling laws: V ∝ r₀² f² (M/Q) q_stable
5. ✅ Practical design methodology
6. ✅ Important negative result: Li+/Na+ infeasible (9 orders of magnitude voltage gap)
7. ✅ Positive result: Microsphere separation feasible (SI = 2.3)

## Timeline

- **Weeks 1-2**: Structure and outline review
- **Weeks 3-4**: Write core theory sections
- **Weeks 5-6**: Create general framework content
- **Weeks 7-8**: Design guidelines and charts
- **Weeks 9-10**: Discussion and polish
- **Weeks 11-12**: Final polish and submit

**Target submission**: 3 months from start

## Notes

- This paper establishes priority on aqueous QIT framework
- Negative result (Li+/Na+ infeasibility) is valuable - prevents wasted research
- Cross-references Paper 2 (microsphere experiments)
- Should be submitted FIRST before Paper 2
