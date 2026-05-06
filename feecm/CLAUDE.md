# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

FEECM (Finite Element Electrochemistry) is a collection of MOOSE-based applications for modeling battery and electrochemical systems. It contains four applications:

- **ec_beta**: Electrochemistry beta version — electrochemistry only (Butler-Volmer, ion transport, current density). Easiest entry point.
- **ecm_test**: Electro-chemo-mechanics test version — adds mechanics coupling on top of electrochemistry.
- **tecm_test**: Thermal-electro-chemo-mechanics — most comprehensive, also adds thermal physics. Active development; forked from `eel` and extended.
- **eel**: Parallel finite-element code for battery physics (mass transport, electrostatics, mechanical deformation, thermal). Originally from literature (UChicago Argonne, 2023); serves as the baseline for `tecm_test`.

## Repository Layout

```
feecm/
├── ec_beta/      ← app source + build artifacts
├── ecm_test/
├── eel/
├── tecm_test/
└── projects/     ← all input files (.i) and project work, separated from app source
    ├── problems-ecBeta/
    ├── problems-ecmTest/
    ├── problems-eel/
    └── problems-tecmTest/
```

The `projects/` folder is the science workspace; the four app folders are the tooling. Each app folder still has its own `doc/` for documenting the app itself; project-level scientific writeups live under `projects/<name>/doc/`.

## Build System

This is a MOOSE framework-based project. Each application has its own Makefile that follows MOOSE conventions:

### Building Applications
```bash
# Build specific applications
cd ec_beta && make -j N
cd ecm_test && make -j N  
cd eel && make -j N

# Where N is the number of parallel jobs
```

### Build Dependencies
- Requires MOOSE framework installation
- Uses standard MOOSE build system with `make`
- Applications link against MOOSE modules (specified in each Makefile)

## Application Structure

Each application follows MOOSE app structure:
- `src/`: Source code organized by object types (kernels, materials, bcs, etc.)
- `include/`: Header files mirroring src structure
- `test/`: Test framework files
- `doc/`: Documentation describing the app itself
- `<app>-opt`: Compiled optimized executable (note: ecm_test's executable is named `ecm-opt`, not `ecm_test-opt`)

Input files (.i) for simulations and test cases live **outside** the app folders, under the top-level `projects/problems-<app>/` trees.

## Testing

Each application has test runners:
```bash
# Run tests for specific applications
./ec_beta/run_tests
./ecm_test/run_tests
./eel/run_tests

# Unit tests
./ec_beta/unit/run_tests
./ecm_test/unit/run_tests
./eel/unit/run_tests
```

## Running Input Files - Important Workflow

**CRITICAL RULE**: Always run input files from their directory location, not from the application root.

### Correct Workflow:
```bash
# Navigate to input file directory first
cd projects/problems-tecmTest/phase2_edl/
# Then run with relative path to the app's executable
../../../tecm_test/tecm_test-opt -i single_ion.i
```

The relative path is `../../../<app>/<exe>`: three `..` to climb out of `projects/problems-<app>/<name>/` to `feecm/`, then into the app folder.

### Per-app executable paths (from a project subdirectory):

| Project tree | Executable invocation |
|---|---|
| `projects/problems-ecBeta/<name>/`  | `../../../ec_beta/ec_beta-opt -i <file>.i` |
| `projects/problems-ecmTest/<name>/` | `../../../ecm_test/ecm-opt -i <file>.i` |
| `projects/problems-eel/<name>/`     | `../../../eel/eel-opt -i <file>.i` |
| `projects/problems-tecmTest/<name>/`| `../../../tecm_test/tecm_test-opt -i <file>.i` |

If a project lives more than one level deep under its `problems-<app>/` tree (e.g., `projects/problems-ecBeta/paper_X/case1/`), add one more `..` for each extra level.

### Why this matters:
- MOOSE resolves relative paths (output directories, file references) relative to current working directory
- Input files often use relative paths like `file_base = rst/output_name`
- Ensures proper file organization and output placement

## MOOSE Modules Used

### ec_beta
- CHEMICAL_REACTIONS: yes
- HEAT_TRANSFER: yes
- PHASE_FIELD: yes
- SOLID_MECHANICS: yes

### ecm_test  
- CHEMICAL_REACTIONS: yes
- CONTACT: yes
- FLUID_PROPERTIES: yes
- HEAT_TRANSFER: yes
- PHASE_FIELD: yes
- POROUS_FLOW: yes
- SOLID_MECHANICS: yes
- XFEM: yes

### eel
- NAVIER_STOKES: yes
- SOLID_MECHANICS: yes

## Key Physics and Components

### Electrochemistry (ec_beta)
- Butler-Volmer kinetics
- Ion transport and migration
- Chemical potential coupling
- Current density calculations

### Electro-Chemo-Mechanics (ecm_test)
- Coupled diffusion-mechanics
- Contact mechanics and delamination
- Stress-assisted transport
- Swelling and growth deformation
- Porous flow modeling

### Battery Modeling (eel)
- Mass transport of charged species
- Electrostatics/Electrodynamics  
- Mechanical deformation
- Thermal effects
- Multi-physics coupling

## Input File Format

Uses MOOSE HIT (Hierarchical Input Text) format with .i extension. Example structure:
- Global parameters and material properties
- Mesh definition
- Variables and equations (kernels)
- Boundary conditions
- Postprocessors and outputs

## Architecture Notes

- Built on MOOSE framework providing FEM discretization
- Massively parallel (MPI + threading)
- Supports adaptive mesh refinement
- Variational framework enables generic coupling implementations
- Object-oriented design with factories for extensibility