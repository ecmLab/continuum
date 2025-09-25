# TECM-no-EN Implementation Guide

This document outlines the implementation of the TECM without electroneutrality (TECM-no-EN) framework based on the derivation in `TECM_noEN_transference.pdf`.

## Overview

The TECM-no-EN framework extends the current EEL-based approach to handle:
- **Explicit valences** (z_i) for each ionic species
- **Space charge effects** via Poisson equation
- **Transference numbers** (t_i) for multi-species transport
- **Electric double layer formation** at interfaces

## Key Equations

### 1. Poisson Equation (Space Charge)
```
-∇ · (ε ∇Ψ) = ρₑ = F(z₊c₊ + z₋c₋) + ρfixed
```

### 2. Species Balance
```
ċᵢ + ∇ · jᵢ = Rᵢ,  i ∈ {+, -}
```

### 3. Species Fluxes
```
jᵢ = -∑ⱼ M⁰ᵢⱼ ∇μⱼ + (tᵢ)/(zᵢF) i
```

### 4. Current Density
```
i = -σ⁰ ∇Φ - (σ⁰/F)(t₊∇μ₊ + t₋∇μ₋)
```

## New Files Created

### 1. Input File
- `examples/paper_ecram/ecram_noEN.i` - ECRAM simulation using TECM-no-EN framework

### 2. C++ Materials
- `include/materials/CurrentDensityNoEN.h/.C` - Multi-species current density calculation
- `include/materials/SpeciesFluxNoEN.h/.C` - Species flux with transference numbers

## Files That Need Modification

### 1. Existing Materials to Update
- `Migration.C` - Update to handle multiple species with transference numbers
- `MassFlux.C` - Extend for multi-species transport
- `ChargeTransferReaction.C` - Update interface conditions for space charge

### 2. New Kernels Needed
- `PoissonEquation.C` - Space charge kernel: `-∇ · (ε ∇Ψ) = -ρₑ`
- `ChargeDensity.C` - Auxiliary kernel: `ρₑ = F(z₊c₊ + z₋c₋)`

### 3. Build System
- Update `CMakeLists.txt` or Makefile to include new materials

## Implementation Steps

### Phase 1: Basic Framework ✓
1. ✅ Create input file template (`ecram_noEN.i`)
2. ✅ Implement `CurrentDensityNoEN` material
3. ✅ Implement `SpeciesFluxNoEN` material
4. ✅ Document implementation approach

### Phase 2: Kernels and AuxKernels
1. ⏳ Create `PoissonEquation` kernel for space charge
2. ⏳ Create `ChargeDensity` auxiliary kernel
3. ⏳ Test basic Poisson equation solution

### Phase 3: Integration
1. ⏳ Update existing materials for multi-species compatibility
2. ⏳ Modify interface kernels for space charge effects
3. ⏳ Test complete ECRAM simulation

### Phase 4: Validation
1. ⏳ Compare with electroneutral limit (ρₑ = 0)
2. ⏳ Verify charge conservation
3. ⏳ Validate Debye length resolution

## Key Implementation Details

### Variable Structure
```cpp
// Primary Variables
c_Li      // Li⁺ concentration [mol/m³]
c_ClO4    // ClO₄⁻ concentration [mol/m³]
Phi       // Transport potential [V]
Psi       // Poisson potential [V]

// Auxiliary Variables
rho_e     // Charge density [C/m³]
E_x, E_y  // Electric field components [V/m]
```

### Material Properties
```cpp
// Transport Properties
sigma_e           // Electrical conductivity [S/m]
D_Li, D_ClO4     // Diffusivities [m²/s]
t_plus, t_minus  // Transference numbers [-]

// Electrochemical Properties
mu_Li, mu_ClO4   // Chemical potentials [J/mol]
i_total          // Current density [A/m²]
j_Li_total       // Li⁺ flux [mol/(m²·s)]
j_ClO4_total     // ClO₄⁻ flux [mol/(m²·s)]
```

### Physical Constants
```cpp
F = 96485         // Faraday constant [C/mol]
R = 8.314         // Gas constant [J/(mol·K)]
z_Li = +1         // Li⁺ valence
z_ClO4 = -1       // ClO₄⁻ valence
epsilon = 8.85e-12 // Permittivity [F/m]
```

## Mesh Requirements

### Debye Length Resolution
The Debye length must be resolved for accurate space charge calculation:
```
λ_D = √(εRT/(F²c)) ≈ 1-10 nm for concentrated electrolytes
```

**Recommendation**: Use mesh refinement near interfaces with element size ≤ λ_D/3.

### Interface Treatment
- Use `BreakMeshByBlockGenerator` to create interfaces
- Apply specialized interface kernels for charge transfer
- Ensure consistent potential definitions across interfaces

## Testing Strategy

### 1. Unit Tests
- Test individual materials with known analytical solutions
- Verify transference number constraints (t₊ + t₋ = 1)
- Check charge conservation (∇ · i = 0 in steady state)

### 2. Benchmark Problems
- 1D planar capacitor with known EDL solution
- Comparison with Gouy-Chapman-Stern theory
- Electroneutral limit recovery (ρₑ → 0)

### 3. ECRAM Validation
- Compare with existing EEL-based results in bulk regions
- Verify EDL formation at electrolyte/electrode interfaces
- Check memory behavior and retention characteristics

## Expected Benefits

### Physical Realism
- ✅ Captures electric double layer formation
- ✅ Includes space charge effects
- ✅ Proper valence handling for multi-valent systems
- ✅ Realistic interface charging

### Numerical Advantages
- ✅ Eliminates electroneutrality constraint singularities
- ✅ More robust convergence for interface problems
- ✅ Natural handling of charge transfer reactions

### Application Scope
- ✅ ECRAM and other memory devices
- ✅ Electrochemical capacitors
- ✅ Battery interfaces with significant EDL effects
- ✅ Solid electrolyte systems

## References

1. `TECM_noEN_transference.pdf` - Mathematical derivation
2. Newman & Thomas-Alyea - "Electrochemical Systems" - Concentrated solution theory
3. Bazant et al. - "Diffuse-charge dynamics in electrochemical systems" - EDL theory