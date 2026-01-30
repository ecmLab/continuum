# ECRAM Device Modeling with TECM Framework

This directory contains the implementation of electrochemical random access memory (ECRAM) device simulation using the TECM (Thermo-Electro-Chemo-Mechanical) framework.

## Physical System

**Device Structure**: 3-terminal side-gated ECRAM with Li+ intercalation into graphene/MoS2 2D heterostructure

**Key Physics**:
- **Electro**: Electric field distribution from gate bias
- **Chemo**: Li+ transport in PEO electrolyte and intercalation into 2D materials
- **Mechanical**: Stress from intercalation-induced swelling

## Mechanism

1. **Gate Bias Applied**: Positive voltage drives Li+ migration toward channel
2. **Li+ Intercalation**: Li+ ions intercalate into van der Waals gaps of 2D heterostructure
3. **Electron Supply**: Electrons provided from Au gate electrode through external circuit
4. **Charge Balance**: ClO4- anions accumulate at gate/electrolyte interface
5. **Memory Effect**: Intercalated Li+ changes channel conductivity (non-volatile)

## Model Implementation

### Core Files
- `ecram_basic.i` - Basic TECM input file for ECRAM simulation
- `README.md` - This documentation

### Key Components Used
- **Variables**: Electric potential (Phi), Li+ concentration (c), displacements (disp_x/y)
- **Materials**: 
  - `BulkChargeTransport` - Electric field distribution
  - `Migration` + `MassDiffusion` - Li+ transport
  - `ChargeTransferReaction` - Intercalation kinetics
  - `SwellingDeformationGradient` - Mechanical swelling
  - `NeoHookeanSolid` - Elastic response

### Physics Coupling
```
Electric Field → Li+ Migration → Intercalation → Swelling → Stress
     ↑                                               ↓
Gate Bias ←------------------- Mechanical Feedback ----
```

## Next Steps

1. **Validate Transport**: Test Li+ diffusion and migration in PEO electrolyte
2. **Calibrate Kinetics**: Fit intercalation exchange current density to experimental data
3. **3D Extension**: Extend to full 3D geometry with realistic device dimensions
4. **Multi-cycle**: Implement write/erase cycling with degradation effects
5. **Experimental Validation**: Compare with device measurements

## Parameters to Tune

### Electrochemical
- `i0` - Exchange current density for intercalation (fit to experiment)
- `U_eq` - Equilibrium potential for Li+ intercalation
- `sigma_e` - Ionic conductivity of PEO electrolyte
- `D_Li` - Li+ diffusivity in PEO

### Mechanical  
- `Omega` - Molar volume of intercalated Li
- `beta` - Swelling coefficient for 2D materials
- `E_2d` - Elastic modulus of graphene/MoS2

### Geometry
- Device dimensions from experimental setup
- Electrolyte thickness and channel length
- 2D material layer thickness

## References

1. User's draft manuscript: "Understanding the Mechanism of Electrochemical Intercalation into Graphene and MoS2 2D Heterostructure"
2. Hu et al. (2023) - EEL theoretical framework 
3. TECM framework implementation in `/Users/howardtu/Documents/modeling/continuum/feecm/tecm_test`