# TECM_TEST Package Structure Analysis

## Overview

**TECM_TEST** (Thermo-Electro-Chemo-Mechanical Test) is a MOOSE-based finite element application for modeling coupled electrochemical and mechanical systems, particularly focused on battery and electrochemical memory device applications.

## Directory Structure

```
tecm_test/
├── .claude/                     # Claude Code configuration
├── .jitcache/                   # JIT compilation cache
├── .libs/                       # Library files
├── build/                       # Build directory
│   ├── header_symlinks/         # Header file symlinks
│   └── unity_src/               # Unity build sources
├── doc/                         # Documentation
│   ├── content/                 # Documentation content
│   └── referencePapers/         # Reference papers
├── examples/                    # Example simulations
│   ├── paper_ecram/             # ECRAM device examples
│   └── tests/                   # Test cases
├── include/                     # Header files
│   ├── auxkernels/             # Auxiliary kernel headers
│   ├── base/                   # Base class headers
│   ├── bcs/                    # Boundary condition headers
│   ├── interfacekernels/       # Interface kernel headers
│   ├── kernels/                # Kernel headers
│   ├── materials/              # Material headers
│   └── utils/                  # Utility headers
├── lib/                        # Library directory
├── rst/                        # ReStructuredText files
├── scripts/                    # Build/utility scripts
├── src/                        # Source code
│   ├── auxkernels/             # Auxiliary kernels
│   ├── base/                   # Base classes
│   ├── bcs/                    # Boundary conditions
│   ├── interfacekernels/       # Interface kernels
│   ├── kernels/                # Primary kernels
│   ├── materials/              # Material models
│   └── main.C                  # Main application file
├── test/                       # Test framework
│   ├── include/                # Test headers
│   ├── lib/                    # Test libraries
│   ├── src/                    # Test source
│   └── tests/                  # Unit tests
├── unit/                       # Unit tests
│   ├── include/                # Unit test headers
│   └── src/                    # Unit test source
├── tecm_test-opt               # Compiled executable
├── tecm_test.yaml              # Configuration file
├── run_tests                   # Test runner script
├── Makefile                    # Build configuration
└── LICENSE                     # License file
```

## Source Code Organization

### Kernels (10 files)
Primary physics equations and governing equations:
- `ADNernstPlanckConvection.C` - Electrochemical ion transport with convection
- `DarcyAdvection.C` - Fluid flow advection using Darcy's law
- `DarcyPressure.C` - Pressure equations for porous media flow
- `EnergyBalanceTimeDerivative.C` - Time derivative for energy conservation
- `MaterialSource.C` - Material source terms
- `NonDimensionalDiffusion.C` - Non-dimensional diffusion equations
- `PoissonEquation.C` - Electrostatic potential equations
- `RankOneDivergence.C` - Vector field divergence operations
- `RankTwoDivergence.C` - Tensor field divergence operations

### Materials (30 files)
Material property definitions and constitutive models:
- **Electrochemical Properties:**
  - `BulkChargeTransport.C` - Bulk ion transport properties
  - `ChargeTransferReaction.C` - Electrode reaction kinetics
  - `ChemicalPotential.C` - Chemical potential calculations
  - `CurrentDensity.C` - Current density calculations
  - `CurrentDensityNoEN.C` - Current density without electro-neutrality

- **Mechanical Properties:**
  - `DeformationGradient.C` - Deformation tensor calculations
  - `SwellingDeformationGradient.C` - Swelling-induced deformation
  - `StressAssisted*.C` - Stress-assisted transport phenomena
  - `ThermalDeformationGradient.C` - Thermal expansion effects

- **Thermodynamic Properties:**
  - `ChemicalEnergyDensity.C` - Chemical energy contributions
  - `ThermalEnergyDensity.C` - Thermal energy calculations
  - `ThermodynamicForce.C` - Driving forces for transport

### Auxiliary Kernels (2 files)
Secondary calculations and derived quantities:
- `ChargeDensity.C` - Charge density calculations
- `DarcyVelocity.C` - Velocity field from Darcy's law

### Boundary Conditions
Interface and boundary condition implementations for coupled physics.

## MOOSE Modules Configuration

From `Makefile`, the application uses:
- **CHEMICAL_REACTIONS**: Chemical reaction modeling
- **HEAT_TRANSFER**: Thermal effects
- **NAVIER_STOKES**: Fluid dynamics
- **SOLID_MECHANICS**: Mechanical deformation
- **THERMAL_HYDRAULICS**: Coupled thermal-hydraulic effects

## Example Applications

### 1. ECRAM Device Modeling (`examples/paper_ecram/`)
Electrochemical Random Access Memory device simulation:
- **Physics**: Li+ intercalation into 2D heterostructures
- **Key Files**: `ecram_full.i`, `ecram_ec.i`
- **Focus**: Gate-controlled Li+ transport and memory effects

### 2. Test Cases (`examples/tests/`)
Various validation and demonstration cases:
- `ex1_constantDisp_2D.i` / `ex1_constantDisp_3D.i` - Constant displacement tests
- `ex2_thermomech.i` - Thermo-mechanical coupling
- `ex3_stress_assisted_diffusion.i` - Stress-assisted transport
- `ex4_electro_chemo_mechanics.i` - Full electro-chemo-mechanics coupling
- `nondimensional_diffusion_demo.i` - Non-dimensional analysis
- `test_thermalmechanics.i` - Thermal mechanics validation

## Key Physics Capabilities

### Electrochemical
- **Ion Transport**: Nernst-Planck equations with migration and diffusion
- **Electrode Kinetics**: Butler-Volmer and charge transfer reactions
- **Electrostatics**: Poisson equation for electric potential
- **Current Density**: Ohmic and Faradaic current calculations

### Mechanical
- **Deformation**: Large deformation finite element framework
- **Swelling**: Volume change from ion intercalation/deintercalation
- **Stress Effects**: Stress-assisted transport and mechanical coupling
- **Thermal Expansion**: Temperature-dependent deformation

### Thermal
- **Heat Transfer**: Conduction and convection
- **Energy Balance**: Thermal energy conservation
- **Coupling**: Temperature effects on transport and mechanical properties

### Multi-Physics Coupling
- **Electro-Chemo**: Electric field effects on ion transport
- **Chemo-Mechanical**: Concentration-dependent swelling/stress
- **Thermo-Mechanical**: Temperature-dependent mechanical properties
- **Full TECM**: All physics coupled simultaneously

## Build and Test Configuration

### Build System
- **Framework**: MOOSE-based with standard Makefile
- **Compiler**: Clang with optimization (`opt` mode)
- **Dependencies**: MOOSE framework installation required

### Testing
- **Test Runner**: `./run_tests` for integration tests
- **Unit Tests**: `./unit/run_tests` for component testing
- **Validation**: Extensive example cases with reference solutions

## Application Focus Areas

1. **Battery Modeling**: Li-ion cell electrochemistry and mechanics
2. **Electrochemical Memory**: ECRAM and other memory device physics
3. **Coupled Transport**: Multi-physics phenomena in electrochemical systems
4. **Material Science**: Stress effects on electrochemical processes
5. **Device Simulation**: Realistic geometry and operating conditions

## Development Status

Based on CLAUDE.md context:
- **Status**: Test version (not yet passed testing phase)
- **Relationship**: Part of FEECM suite alongside `ec_beta` and `eel`
- **Focus**: Electro-chemo-mechanics coupling with advanced features
- **Applications**: Research-grade simulations for battery and memory devices

## Dependencies and Requirements

- **MOOSE Framework**: Primary computational framework
- **MPI**: Parallel execution support
- **PETSc**: Linear/nonlinear solver backend
- **libMesh**: Finite element discretization
- **Compiler**: Modern C++ with C++14/17 support

This structure analysis reveals TECM_TEST as a comprehensive multi-physics simulation platform specialized for electrochemical systems with strong mechanical coupling capabilities.