# Two-Interface Double Layer Model Proposal

## Model Overview

**Geometry:** 1D interval with TWO electrode/electrolyte interfaces

```
[Electrode L] |<-- DDL -->| [Electrolyte] |<-- DDL -->| [Electrode R]
    x = 0              x = L/2              x = L
  Boundary 1         (center)           Boundary 2
```

## Configuration Options

### **Option 1: Symmetric Capacitor (Recommended)**
- **Left electrode (x=0)**: Potential = +V/2
- **Right electrode (x=L)**: Potential = -V/2
- **Center (x=L/2)**: Symmetry plane, E-field = 0
- **Physics**: Two identical double layers forming back-to-back

**Advantages:**
- Physically realistic (like a capacitor)
- Symmetric solution reduces numerical issues
- Can use symmetry boundary condition
- Clear interpretation of results

### **Option 2: Asymmetric Configuration**
- **Left electrode**: Potential = +V₁
- **Right electrode**: Potential = +V₂ (V₁ ≠ V₂)
- **Physics**: Two different double layers

**Advantages:**
- More general
- Can study asymmetric effects
- Useful for battery-like systems

### **Option 3: One Grounded, One Biased**
- **Left electrode**: Potential = V
- **Right electrode**: Ground (0 V)
- **Physics**: Similar to single interface but with two boundaries

## Recommended Configuration: Option 1 (Symmetric)

### Geometry
```matlab
L_cell = xD * 200  % Total length: 200 Debye lengths
% Electrode positions: x = 0 and x = L_cell
```

### Boundary Conditions

**At x = 0 (Left Electrode):**
1. **Electrostatics:**
   - Surface charge density: `rho_surf_L = ε₀*ε_S*(phiM_L - phi)/xS`
   - phiM_L = +phiM/2 = +5 mV

2. **Transport:**
   - Zero flux (blocking electrode)
   - No faradaic reactions

**At x = L_cell (Right Electrode):**
1. **Electrostatics:**
   - Surface charge density: `rho_surf_R = ε₀*ε_S*(phiM_R - phi)/xS`
   - phiM_R = -phiM/2 = -5 mV

2. **Transport:**
   - Zero flux (blocking electrode)

**In Domain:**
- Coupled Poisson-Nernst-Planck equations
- Initial condition: uniform bulk concentration

### Parameters

Same as one-interface model:
- `DA = 1e-12 m²/s`
- `DX = 2*DA`
- `cA_bulk = 10 mol/m³`
- `eps_PEO = 10`
- `phiM = 10 mV` (total potential difference)
- `T0 = 25°C`

**Modified:**
- `L_cell = xD * 200` (longer to accommodate two double layers)
- `phiM_L = phiM/2 = +5 mV`
- `phiM_R = -phiM/2 = -5 mV`

### Expected Results

1. **Spatial profiles at steady state:**
   - Electric potential: S-shaped, drops from +5 mV to -5 mV
   - Cation concentration: depleted near x=0, accumulated near x=L
   - Anion concentration: accumulated near x=0, depleted near x=L
   - Two distinct double layers visible

2. **Time evolution:**
   - Initially: uniform concentrations
   - Transient: Double layers form simultaneously at both ends
   - Steady state: Symmetric concentration/potential profiles
   - Current: Decays to zero (capacitor fully charged)

3. **Key physics:**
   - Total capacitance = C_L series with C_R
   - For symmetric case: C_total = C_single/2
   - Charge conservation: Q_L + Q_R = 0

## Implementation Changes from One-Interface Model

### 1. Geometry
- Keep interval but interpret differently
- x = 0: left electrode
- x = L: right electrode

### 2. Parameters
```matlab
model.param.set('phiM_L', '+phiM/2', 'Left electrode potential vs PZC');
model.param.set('phiM_R', '-phiM/2', 'Right electrode potential vs PZC');
model.param.set('L_cell', 'xD*200', 'Total cell length (longer for 2 interfaces)');
```

### 3. Variables
```matlab
% Left electrode
var1.set('deltaphi_L', 'phiM_L-phi', 'Left electrode-OHP potential difference');
var1.set('rho_surf_L', 'epsilon0_const*eps_S*deltaphi_L/xS', 'Left surface charge density');

% Right electrode
var1.set('deltaphi_R', 'phiM_R-phi', 'Right electrode-OHP potential difference');
var1.set('rho_surf_R', 'epsilon0_const*eps_S*deltaphi_R/xS', 'Right surface charge density');
```

### 4. Electrostatics BCs
```matlab
% Left electrode (boundary 1)
model.component('comp1').physics('es').create('sfcd1', 'SurfaceChargeDensity', 0);
model.component('comp1').physics('es').feature('sfcd1').selection.set([1]);
model.component('comp1').physics('es').feature('sfcd1').set('rhoqs', 'rho_surf_L');

% Right electrode (boundary 2)
model.component('comp1').physics('es').create('sfcd2', 'SurfaceChargeDensity', 0);
model.component('comp1').physics('es').feature('sfcd2').selection.set([2]);
model.component('comp1').physics('es').feature('sfcd2').set('rhoqs', 'rho_surf_R');

% Remove ground condition - both electrodes have surface charge BCs
```

### 5. Transport BCs
Both electrodes already have zero flux (no concentration BC needed)
- Default no-flux condition at boundaries applies

### 6. Mesh
Refine near BOTH electrodes:
```matlab
% Refinement at x=0 (boundary 1)
mesh.feature('edg1').feature('size_left').selection.set([1]);
mesh.feature('edg1').feature('size_left').set('hmax', 'h_max_surf');

% Refinement at x=L (boundary 2)
mesh.feature('edg1').feature('size_right').selection.set([2]);
mesh.feature('edg1').feature('size_right').set('hmax', 'h_max_surf');
```

### 7. Data Extraction
Extract at BOTH electrodes:
- Left electrode: boundary 1 (x=0)
- Right electrode: boundary 2 (x=L)

## Questions for You

1. **Which configuration do you prefer?**
   - Option 1 (Symmetric: +V/2 and -V/2)?
   - Option 2 (Asymmetric)?
   - Option 3 (One grounded)?

2. **Domain length:**
   - 200×xD (my suggestion for clear separation)?
   - Different factor?

3. **Additional outputs:**
   - Extract data at center point (x=L/2)?
   - Compare left vs right double layer properties?
   - Plot total charge Q_L + Q_R to verify conservation?

Please confirm your preference and I'll generate the complete ddl_twoInterface.m file!
