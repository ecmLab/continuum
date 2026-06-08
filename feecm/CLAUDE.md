# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

FEECM (Finite Element Electrochemistry) is a collection of MOOSE-based applications for modeling battery and electrochemical systems. It contains two applications:

- **ecm_test**: Electro-chemo-mechanics — covers electrochemistry (Butler-Volmer, ion transport, current density), mechanics coupling (diffusion-mechanics, contact, swelling, porous flow), and consolidates the former `ec_beta` capabilities.
- **tecm_test**: Thermal-electro-chemo-mechanics — battery physics (mass transport, electrostatics, mechanical deformation, thermal) plus phase-oriented EDL development and a non-dimensional layer. Originally forked from `eel` (UChicago Argonne literature baseline, 2023); now consolidates the full `eel` capability.

> **Historical notes:**
> - The repo previously had an `ec_beta` app (electrochemistry-only beta version). It was consolidated into `ecm_test` and removed; all 74 of its registered classes now live in `ecm_test` and the `ecm-opt` executable is a strict superset of the former `ec_beta-opt`.
> - The repo previously had a standalone `eel` app (parallel FE battery physics). It was consolidated into `tecm_test` and removed; all of its registered classes now live in `tecm_test` and the `tecm_test-opt` executable is a strict superset of the former `eel-opt`. The former `projects/eel/` work moved under `projects/tecmTest/`.

## Repository Layout

```
feecm/
├── ecm_test/
├── tecm_test/
├── patchMoose2Feecm/  ← MOOSE framework patches (manual application)
└── projects/          ← all input files (.i) and project work, separated from app source
    ├── ecmTest/        (← merged from former ecBeta + original ecmTest)
    └── tecmTest/       (← original tecmTest + merged former projects/eel)
```

The `projects/` folder is the science workspace; the two app folders are the tooling. Each app folder still has its own `doc/` for documenting the app itself; project-level scientific writeups live under `projects/<name>/doc/`.

## Build System

This is a MOOSE framework-based project. Each application has its own Makefile that follows MOOSE conventions:

### Building Applications
```bash
# Build specific applications
cd ecm_test && make -j N
cd tecm_test && make -j N

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

Input files (.i) for simulations and test cases live **outside** the app folders, under the top-level `projects/<app>/` trees.

## Testing

Each application has test runners:
```bash
# Run tests for specific applications
./ecm_test/run_tests
./tecm_test/run_tests

# Unit tests
./ecm_test/unit/run_tests
./tecm_test/unit/run_tests
```

## Running Input Files - Important Workflow

**CRITICAL RULE**: Always run input files from their directory location, not from the application root.

### Correct Workflow:
```bash
# Navigate to input file directory first
cd projects/tecmTest/phase2_edl/
# Then run with relative path to the app's executable
../../../tecm_test/tecm_test-opt -i single_ion.i
```

The relative path is `../../../<app>/<exe>`: three `..` to climb out of `projects/<app>/<name>/` to `feecm/`, then into the app folder.

### Per-app executable paths (from a project subdirectory):

| Project tree | Executable invocation |
|---|---|
| `projects/ecmTest/<name>/`  | `../../../ecm_test/ecm-opt -i <file>.i` |
| `projects/tecmTest/<name>/` | `../../../tecm_test/tecm_test-opt -i <file>.i` |

If a project lives more than one level deep under its `<app>/` tree (e.g., `projects/ecmTest/paper_fengyu_anodeLess_BL/randomPores/`), add one more `..` for each extra level.

### Why this matters:
- MOOSE resolves relative paths (output directories, file references) relative to current working directory
- Input files often use relative paths like `file_base = rst/output_name`
- Ensures proper file organization and output placement

## MOOSE Modules Used

### ecm_test
- CHEMICAL_REACTIONS: yes
- CONTACT: yes
- FLUID_PROPERTIES: yes
- HEAT_TRANSFER: yes
- PHASE_FIELD: yes
- POROUS_FLOW: yes
- SOLID_MECHANICS: yes
- XFEM: yes

### tecm_test
- CHEMICAL_REACTIONS: yes
- HEAT_TRANSFER: yes
- NAVIER_STOKES: yes
- SOLID_MECHANICS: yes
- THERMAL_HYDRAULICS: yes
- (see `tecm_test/Makefile` for the authoritative list)

## Key Physics and Components

### Electro-Chemo-Mechanics (ecm_test)
- Butler-Volmer kinetics, ion transport, migration, current density (former ec_beta capabilities)
- Coupled diffusion-mechanics
- Contact mechanics and delamination
- Stress-assisted transport
- Swelling and growth deformation
- Porous flow modeling
- MIEC (mixed ionic/electronic conductors), Nernst-Planck convection, anisotropic diffusion / conductivity

### Thermal-Electro-Chemo-Mechanics (tecm_test)
- Battery physics: mass transport of charged species, electrostatics/electrodynamics, mechanical deformation, thermal effects, multi-physics coupling (the former `eel` capability)
- Phase-oriented EDL development and non-dimensional layer
- See `tecm_test/AGENTS.md` for the phase-oriented development conventions

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
