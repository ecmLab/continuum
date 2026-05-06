# ddl1.m Update Notes

## Summary
Updated `ddl1.m` to match parameters from `ddl0.m` and added time-dependent simulation capability.

## Key Changes

### 1. **Parameters Updated to Match ddl0.m**

| Parameter | Old Value | New Value | Notes |
|-----------|-----------|-----------|-------|
| `DA` | `1e-9[m^2/s]` | `1e-12[m^2/s]` | Diffusion coefficient (1000x slower - solid electrolyte) |
| `DX` | `DA` | `DA/5` | Anion diffuses 5x slower than cation |
| `cA_bulk` | `10[mol/m^3]` | `110[mol/m^3]` | 11x higher bulk concentration |
| `cX_bulk` | `cA_bulk` | `cA_bulk` | (same formula, but value changes) |
| `eps_PEO` | N/A | `5.0` | Added: Relative permittivity of PEO electrolyte |
| `phiM` | `10[mV]` | `1[mV]` | Back to 1 mV electrode potential |

### 2. **Study Configuration**

**Old:**
- Single stationary study only

**New:**
- Stationary study (inactive, for parametric sweeps if needed)
- **Time-dependent study (active)**
  - Time range: 0.001 to 10 ms
  - Logarithmic time stepping: `10^{range(log10(0.001),1/3,log10(10))}`
  - Time unit: milliseconds

### 3. **Mesh Configuration**

**Old:**
- Automatic mesh with `autoMeshSize(5)`

**New:**
- Manual mesh matching ddl0.m exactly:
  - Edge meshing with two size controls
  - `size1`: General mesh with `hmax = L_cell/20`
  - `size2`: Refined mesh at electrode (point 1) with `hmax = xD/100`

### 4. **Data Extraction**

**Old:**
- Only spatial profiles at final state

**New:**
- **Spatial profiles** at final time step (x-direction):
  - `x_data`, `phi_data`, `cA_data`, `cX_data`, `Es_data`, `rho_data`

- **Time evolution** at electrode position (x=0):
  - `t_data`: time points (ms)
  - `phi_t`: potential vs time at electrode
  - `cA_t`: cation concentration vs time
  - `cX_t`: anion concentration vs time

### 5. **Output Files**

**Old:**
- `ddl_data.csv` - single CSV with spatial data

**New:**
- `ddl_data_spatial.csv` - Spatial profiles at final time
- `ddl_data_time.csv` - Time evolution at electrode
- Both saved in addition to `.mat` file

### 6. **Plots Generated**

**Old (4 plots):**
1. Electric potential vs distance
2. Ion concentrations vs distance
3. Electric field vs distance
4. Space charge density vs distance

**New (5 plots):**
1. Electric potential vs distance (final state)
2. Ion concentrations vs distance (final state)
3. Electric field vs distance (final state)
4. Space charge density vs distance (final state)
5. **NEW: Time evolution at electrode** (3 subplots):
   - Potential vs time (semilog)
   - Cation concentration vs time (semilog)
   - Anion concentration vs time (semilog)

## Physical Interpretation

### Changed Parameters Context

**Slower diffusion (DA = 1e-12 m²/s):**
- Typical for solid polymer electrolytes (PEO)
- ~1000x slower than aqueous systems
- Justifies time-dependent study (longer equilibration times)

**Higher concentration (110 mol/m³):**
- Typical for concentrated battery electrolytes
- Shorter Debye length (~0.9 nm vs ~3 nm)
- Stronger screening effects

**Asymmetric diffusion (DX = DA/5):**
- Reflects different ion mobilities
- Anions typically larger/slower in PEO
- Creates interesting time-dependent concentration asymmetry

### Time-Dependent Physics

The time-dependent study captures:
1. **Initial transient:** How double layer forms from initial uniform state
2. **Charge relaxation:** Time to establish equilibrium potential
3. **Diffusion timescales:** Different rates for cations vs anions
4. **Steady state:** Final equilibrium matching analytical Gouy-Chapman theory

Characteristic time: τ ~ xD²/D ~ (1e-9)²/(1e-12) ~ 1 ms

Time range (0.001-10 ms) spans from early transient to full equilibration.

## Running the Updated Script

```matlab
% In MATLAB with COMSOL Server running:
cd /Users/howardtu/Documents/modeling/continuum/feecm/comsol/ddl
ddl1

% Expected runtime: ~30-60 seconds (time-dependent solve is slower)
```

## Validation

To validate results:
1. Check that steady state (t=10 ms) matches Gouy-Chapman theory
2. Verify exponential decay with Debye length ~0.9 nm
3. Confirm anion/cation asymmetry due to different diffusivities
4. Compare with ddl0.mph results if available

## Files Modified

- `ddl1.m` - Main simulation script (updated)

## Files Created

- `UPDATE_NOTES.md` - This file
- `validate_results.m` - Result validation script (from previous session)

## 2025-10-09 — Phase 1 Current Extraction Fix

- In `oneInterface.m`, switched the current definitions to use the built-in normal total flux variables (`tds.ntflux_cA`, `tds.ntflux_cX`). This measures the ionic current crossing the electrode surface with the correct orientation (positive toward +x).
- Regenerate the Phase 1 outputs after pulling this change to update the CSVs/plots with the corrected current directions.

## 2025-10-09 — Phase 4 EIS Implementation Complete

### sinusoidal.m Enhancement
- Merged `analyze_impedance.m` coherently into `sinusoidal.m`
- Integrated FFT-based impedance analysis (lines 345-466)
- Consolidated all plotting into single section
- Removed redundant figures
- Generates 6 comprehensive figures including Nyquist and Bode plots

### eis.m - New Multi-Frequency EIS Implementation
- Created based on working `sinusoidal.m` code (not old broken EIS code)
- **Key achievement**: Successfully captures Warburg diffusion impedance (45° line in Nyquist plot)

**Critical parameters for Warburg visibility:**
- Reduced diffusion coefficients: `DA = 1e-14 m²/s` (100× slower than realistic)
- Frequency range: 1000 Hz to ~4 Hz (15 points, logarithmic)
- Excitation: 10 mV AC, 0 V DC
- Simulation: 10 cycles per frequency, 50 points per cycle

**Implementation highlights:**
- Midpoint current extraction (avoids Stern layer complications)
- FFT-based impedance with Hanning window
- Frequency sweep from high to low (standard EIS convention)
- Bode plots with reversed frequency axis (high freq on left)

**Results:**
- Solution resistance: Rs ≈ 0.192 Ω·m² (vs 0.135 Ω·m² theoretical, 40% higher due to non-uniform conductivity)
- Clear 45° Warburg line in Nyquist plot at mid-frequencies (10-100 Hz)
- Transition from resistive (high freq) → diffusion (mid freq) → capacitive (low freq)

**Output files:**
- `rst/eis/eis_data.csv` - Frequency, Z', Z", |Z|, phase
- `rst/eis/eis_results.mat` - Full MATLAB workspace
- `rst/eis/eis_nyquist.png` - Nyquist plot showing Warburg line
- `rst/eis/eis_bode.png` - Bode magnitude and phase plots
- `rst/eis/eis_solved.mph` - COMSOL model with solution

### Documentation Updates
- `doc/EIS_ISSUES_AND_FIXES.md` - Comprehensive development history, physics explanation, lessons learned
- `doc/EIS_Documentation.tex` - Updated LaTeX documentation with current implementation details

### Key Lessons Learned
1. Don't debug broken code - reuse working code (eis.m based on sinusoidal.m)
2. Midpoint current extraction avoids boundary complications
3. Physics must match frequency range (reduced D to make Warburg visible)
4. FFT is more robust than least-squares fitting
5. Simple frequency loop beats complicated adaptive logic

## 2025-10-20 — Phase 5 Reference Paper Implementation (eisRef)

### eis_ref.m - Reference Implementation Following Mei et al. (2018)
- Created new implementation following methodology from "Physical Interpretations of Nyquist Plots for EDLC Electrodes and Devices" (J. Phys. Chem. C 2018, 122, 194-206)
- **Key difference from eis.m**: Uses paper's exact methodology for validation and comparison

**Implementation details:**
- **Single electrode configuration** (blocking electrode at x=0, bulk electrolyte at x=L)
- **Boundary conditions** following paper eq. S.5a, S.7-S.10:
  - x=0: Surface charge density BC (Stern layer), no flux for species
  - x=L: Ground potential, bulk concentration
- **Least-squares fitting** for impedance extraction (paper method, not FFT):
  - Fit I(t) = I_dc + I_sin·sin(ωt) + I_cos·cos(ωt)
  - Z = V_ac / (I_sin - j·I_cos)
  - Track R² goodness-of-fit for quality control
- **Broader frequency range**: 0.1 Hz to 1 MHz (20 points) vs eis.m's 4-1000 Hz
- **Current extraction at electrode boundary** (x=0) instead of midpoint

**Physical configuration:**
- Domain length: L = 100 λ_D (single electrode, vs 200 λ_D for two electrodes in eis.m)
- Same material parameters: D_A = 1e-14 m²/s, D_X = 2e-14 m²/s, c_bulk = 10 mol/m³
- Excitation: 10 mV AC, 0 V DC bias
- 10 cycles per frequency, 50 time points per cycle

**Resistance identification following paper:**
- R_A: Real impedance at highest frequency (electrode resistance)
- R_AB: Semicircle diameter (bulk electrolyte resistance R_∞)
- R_BC: Additional resistance from Warburg/diffusion (R_D)
- Theoretical validation: κ = (F²/RT) Σ z_i² D_i c_i, R_s = L/κ

**Output files:**
- `rst/eisRef/eis_ref_data.csv` - Full results including fit quality (R²)
- `rst/eisRef/eis_ref_results.mat` - MATLAB workspace with all data
- `rst/eisRef/eis_ref_nyquist.png` - Nyquist plot with frequency labels
- `rst/eisRef/eis_ref_bode.png` - Bode magnitude and phase
- `rst/eisRef/eis_ref_fit_quality.png` - Least-squares fit R² vs frequency
- `rst/eisRef/eis_ref_solved.mph` - COMSOL model

**Comparison with eis.m:**
| Feature | eis.m | eis_ref.m |
|---------|-------|-----------|
| Configuration | Two electrodes (symmetric) | Single electrode (blocking) |
| Current extraction | Midpoint (L/2) | Electrode boundary (x=0) |
| Impedance method | FFT with Hanning window | Least-squares fitting |
| Frequency range | 4-1000 Hz (15 pts) | 0.1 Hz-1 MHz (20 pts) |
| Domain length | 200 λ_D | 100 λ_D |
| Focus | Warburg visibility | Full spectrum + validation |
| Validation | Theoretical R_s | Paper methodology (R_A, R_AB, R_BC) |

**Purpose:**
- Validate our MPNP implementation against published paper methodology
- Enable direct comparison with paper's Figure 2a (single electrode)
- Foundation for systematic parameter studies (paper Table 2)
- Explore full frequency spectrum from capacitive to resistive regimes
