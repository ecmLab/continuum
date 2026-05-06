# Phase 5: eisRef - Reference Paper Implementation

## Overview

`eis_ref.m` implements Electrochemical Impedance Spectroscopy (EIS) following the exact methodology from:

> Mei, B.-A., Munteshari, O., Lau, J., Dunn, B., & Pilon, L. (2018).
> **Physical Interpretations of Nyquist Plots for EDLC Electrodes and Devices**
> *J. Phys. Chem. C*, 122, 194-206.
> DOI: [10.1021/acs.jpcc.7b10582](https://doi.org/10.1021/acs.jpcc.7b10582)

This implementation serves as a **validation reference** and **foundation** for systematic parameter studies comparing our MPNP model with published results.

## Key Differences from eis.m

Our original `eis.m` was optimized to visualize Warburg diffusion impedance. This reference implementation (`eis_ref.m`) follows the paper's methodology exactly for validation purposes:

| Feature | eis.m (Warburg-focused) | eis_ref.m (Paper reference) |
|---------|------------------------|------------------------------|
| **Configuration** | Two electrodes (asymmetric) | Single electrode (blocking) |
| **Domain** | 200 λ_D (symmetric cell) | 100 λ_D (semi-infinite) |
| **Current extraction** | Midpoint (x = L/2) | Electrode boundary (x = 0) |
| **Impedance method** | FFT with Hanning window | Least-squares fitting |
| **Frequency range** | 4-1000 Hz (15 points) | 0.1 Hz - 1 MHz (20 points) |
| **Focus** | Warburg line (45°) | Full spectrum + resistance ID |
| **Validation** | Theoretical R_s | R_A, R_AB, R_BC methodology |
| **Output** | Nyquist + Bode | + Fit quality plot |

## Governing Equations (Paper SI)

### 1. Modified Poisson-Nernst-Planck (MPNP) Model

**Poisson equation** (S.2a):
```
∇·(ε₀εᵣ∇φ) = -F Σᵢ zᵢcᵢ
```

**Nernst-Planck equation** (S.2b):
```
∂cᵢ/∂t = -∇·Nᵢ
```

**Ionic flux** (S.3):
```
Nᵢ = -Dᵢ∇cᵢ - (zᵢF/RT)Dᵢcᵢ∇φ
```

### 2. Boundary Conditions (Single Electrode)

**At electrode (x = 0)**:
- **Electric potential** (S.7): Surface charge density
  ```
  σ_Stern = (ε₀εₛ/δₛ)(φₘ - φ|ₓ₌₀)
  ```
- **Species** (S.9): No flux (blocking electrode)
  ```
  Nᵢ·n = 0
  ```

**At bulk (x = L)**:
- **Electric potential** (S.8): Ground
  ```
  φ = 0
  ```
- **Species** (S.10): Bulk concentration
  ```
  cᵢ = c₀,ᵢ
  ```

### 3. Sinusoidal Excitation (S.5a)

```
φₘ(t) = V_DC + V_AC sin(ωt)
```

### 4. Impedance Calculation

**Least-squares fit** of current:
```
I(t) = I_DC + I_sin·sin(ωt) + I_cos·cos(ωt)
```

**Complex impedance**:
```
Z(ω) = V_AC / (I_sin - j·I_cos)
```

## Physical Parameters

### Material Properties
```matlab
T       = 25°C (298.15 K)      % Temperature
D_A     = 1e-14 m²/s           % Cation diffusion coefficient
D_X     = 2e-14 m²/s           % Anion diffusion coefficient
c_bulk  = 10 mol/m³            % Bulk concentration (both ions)
z_A     = +1                   % Cation charge
z_X     = -1                   % Anion charge
ε_PEO   = 10.0                 % Electrolyte permittivity
ε_S     = 10.0                 % Stern layer permittivity
δ_S     = 0.2 nm               % Stern layer thickness
```

### Derived Quantities
```matlab
V_therm = RT/F ≈ 25.7 mV       % Thermal voltage
λ_D     ≈ 2.41 nm              % Debye length
L       = 100 λ_D ≈ 241 nm     % Domain length
```

### Excitation
```matlab
V_DC    = 0.0 V                % DC bias
V_AC    = 0.01 V (10 mV)       % AC amplitude
f       = 0.1 Hz to 1 MHz      % Frequency range (20 points, log-spaced)
```

## Resistance Identification

Following paper Section 3.2 and Figure 3:

### Nyquist Plot Features

```
    -Z" (Ω·m²)
        ↑
        │     ╱╲  ← Semicircle (bulk electrolyte)
        │    ╱  ╲
        │   ╱    ╲
        │  ╱      ╲___
        │ ╱           ╲___
        │╱                ╲___  ← Warburg line (45°)
        A────B───────────────C────→ Z' (Ω·m²)
```

### Key Points
- **Point A** (high frequency): `Z'_A = R_e` (electrode resistance)
- **Point B** (semicircle end): `Z'_B = R_e + R_∞`
- **Point C** (low frequency): `Z'_C = R_e + R_∞ + R_D`

### Derived Resistances
```matlab
R_AB = R_∞     % Bulk electrolyte resistance (semicircle diameter)
R_BC = R_D     % Diffuse layer resistance (Warburg region)
R_B  = R_e + R_∞  % Internal resistance (galvanostatic IR drop)
```

### Theoretical Validation

**Solution conductivity** (S.21):
```
κ = (F²/RT) Σᵢ zᵢ²Dᵢcᵢ
```

**Solution resistance**:
```
R_s = L/κ
```

For our parameters:
```matlab
κ = (96485²/(8.314×298.15)) × [(1²×1e-14 + 1²×2e-14) × 10]
  ≈ 1.166 S/m

R_s = 241e-9 / 1.166 ≈ 2.07e-7 Ω·m²
```

Compare with `R_AB` from simulation to validate implementation.

## Running the Simulation

### Prerequisites
- MATLAB with COMSOL Livelink
- COMSOL Multiphysics 5.x or later

### Execution
```matlab
% Navigate to ddl directory
cd /Users/howardtu/Documents/modeling/continuum/feecm/comsol/ddl

% Run simulation
eis_ref

% Expected runtime: ~5-15 minutes (20 frequencies × ~30s each)
```

### Progress Monitoring
The script provides detailed progress information:
```
=== EIS REFERENCE IMPLEMENTATION ===
Based on: Mei et al., J. Phys. Chem. C 2018, 122, 194-206

EIS Configuration:
  Frequency range: 0.1 Hz to 1.0e+06 Hz (20 points, log-spaced)
  DC bias: 0.000 V
  AC amplitude: 0.010 V (10.0 mV)
  ...

========================================
Frequency 1/20: 1.00e+06 Hz
========================================
  Period: 1.0000e-06 s
  Total time: 1.0000e-05 s (10 cycles)
  ...

  Least-squares fit:
    I_dc = ...
    I_sin = ...
    I_cos = ...
    R^2 = 0.999xxx

  Impedance:
    Z = ... Ohm*m^2
    |Z| = ...
    Phase = ... deg
```

## Output Files

All results are saved in `rst/eisRef/`:

### 1. Data Files
- **`eis_ref_data.csv`**: Complete results table
  ```
  freq_Hz | Z_real_Ohm_m2 | Z_imag_Ohm_m2 | Z_mag_Ohm_m2 | phase_deg |
          | I_dc_A_m2 | I_sin_A_m2 | I_cos_A_m2 | fit_R2
  ```

- **`eis_ref_results.mat`**: MATLAB workspace
  - All arrays: `frequencies`, `Z_real_all`, `Z_imag_all`, etc.
  - Parameters: `V_dc`, `V_ac`, `num_cycles`, `pts_per_cycle`
  - Fit quality: `fit_R2_all`

- **`eis_ref_solved.mph`**: COMSOL model with last frequency solution

### 2. Plots

- **`eis_ref_nyquist.png`**: Nyquist plot
  - Blue line with markers: EIS data
  - Red square: High frequency point (1 MHz)
  - Green square: Low frequency point (0.1 Hz)
  - Frequency labels at 5 key points
  - Equal axis scaling for correct geometry

- **`eis_ref_bode.png`**: Bode plots (2 subplots)
  - Top: |Z| vs frequency (log-log)
  - Bottom: Phase vs frequency (semilog)
  - Both with reversed x-axis (high freq on left)

- **`eis_ref_fit_quality.png`**: Goodness-of-fit
  - R² coefficient vs frequency
  - Quality check for least-squares fitting
  - Should be > 0.99 for all frequencies

### 3. Console Output

Final summary includes:
```
=== RESISTANCE IDENTIFICATION ===
Following Mei et al. (2018) methodology:

R_A (Point A, electrode resistance):
  R_A = ... Ohm*m^2 (at 1e+06 Hz)

R_B (Point B, end of semicircle):
  R_B = ... Ohm*m^2 (at ... Hz)

R_C (Point C, low frequency limit):
  R_C = ... Ohm*m^2 (at 1e-01 Hz)

Derived Resistances:
  R_AB = R∞ (bulk electrolyte): ... Ohm*m^2
  R_BC = R_D (diffuse layer): ... Ohm*m^2
  R_∞ (internal resistance): ... Ohm*m^2

Theoretical Solution Resistance:
  κ (conductivity) = ... S/m
  R_s (theory) = L/κ = ... Ohm*m^2
  R_AB (simulation) = ... Ohm*m^2
  Ratio (sim/theory) = ...
```

## Validation Criteria

### 1. Fit Quality
- R² > 0.99 for all frequencies
- If R² < 0.99, check:
  - Simulation reached steady state (last 2 cycles)
  - Sufficient time points per cycle (≥50)
  - Excitation amplitude not too small

### 2. Nyquist Plot Shape
- **High frequency**: Approaches R_A (electrode resistance)
- **Mid frequency**: Semicircle (bulk electrolyte)
- **Low frequency**: Warburg line or vertical capacitive line

### 3. Theoretical Comparison
- R_AB should match R_s = L/κ within ~30-50%
- Discrepancy expected due to:
  - Non-uniform conductivity (concentration gradients)
  - Double layer effects
  - Stern layer contribution

### 4. Phase Behavior
- **High frequency**: Phase → 0° (resistive)
- **Mid frequency**: Phase between 0° and 90° (mixed)
- **Low frequency**: Phase → 90° (capacitive) or 45° (Warburg)

## Comparison with Paper Results

### Paper Figure 2a (Single Electrode)
Our setup matches paper's single electrode case:
- Same MPNP model
- Same boundary conditions
- Similar parameter range

**Expected features**:
1. Small R_A at high frequency (electrode resistance)
2. Semicircle in Nyquist plot (R_∞)
3. Transition to Warburg or capacitive at low frequency

### Paper Table 2 (Parameter Studies)
Future extensions can vary parameters systematically:
- Diffusion coefficients: D = 10⁻⁹ to 10⁻¹⁵ m²/s
- Concentrations: c₀ = 0.1 to 1000 mol/m³
- Potentials: φ_m = 10 mV to 1 V
- Electrode thickness: L = 10 to 100 nm

## Troubleshooting

### Problem: Fit quality R² < 0.99
**Solution**:
- Increase `num_cycles` (e.g., 15 or 20)
- Increase `pts_per_cycle` (e.g., 100)
- Check convergence at low frequencies

### Problem: Simulation too slow
**Solution**:
- Reduce `num_freqs` (e.g., 15 instead of 20)
- Reduce frequency range (e.g., 1 Hz to 100 kHz)
- Coarsen mesh (increase `h_max`)

### Problem: Strange Nyquist plot shape
**Solution**:
- Check parameter values (D, c, ε)
- Verify boundary conditions
- Ensure L is large enough (≥ 100 λ_D)
- Check for numerical instabilities

### Problem: No Warburg impedance visible
**Expected**: With D = 1e-14 m²/s, Warburg should appear around 10-1000 Hz
**Solutions**:
- Reduce D further (e.g., 1e-15 m²/s)
- Add more frequencies in mid-range (10-1000 Hz)
- Check domain length (L should be comparable to diffusion length)

## Future Directions

### 1. Parametric Studies
Follow paper Table 2 to explore:
- Effect of diffusion coefficient on Warburg frequency
- Effect of concentration on Debye length
- Effect of potential on double layer capacitance

### 2. Two-Electrode System
Extend to symmetric two-electrode cell (paper Figure 2b):
- Different resistance identification
- Device-level impedance
- Comparison with oneInterface.m / twoInterface.m

### 3. Porous Electrodes
Incorporate porosity models:
- Transmission line model
- Distributed resistance and capacitance
- Compare with experimental EDLCs

### 4. Nonlinear Effects
Explore large-amplitude perturbations:
- Harmonic distortion
- Modified Stern layer model
- Ion steric effects

## References

1. **Primary Paper**:
   Mei, B.-A., et al. (2018). "Physical Interpretations of Nyquist Plots for EDLC Electrodes and Devices." *J. Phys. Chem. C*, 122, 194-206.

2. **Supporting Information**:
   Detailed equations and boundary conditions (S.1-S.21)

3. **Related Work**:
   - Gouy-Chapman-Stern theory of double layers
   - Warburg impedance in diffusion-limited systems
   - EIS analysis methods and interpretation

## Contact

For questions about this implementation:
- Check `doc/UPDATE_NOTES.md` for recent changes
- Compare with `eis.m` for Warburg-focused version
- Refer to paper SI for detailed equations

---

**Created**: 2025-10-20
**Status**: Initial implementation
**Next Steps**: Run simulation, validate against paper Figure 2a
