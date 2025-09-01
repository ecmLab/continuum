# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

FEECM (Finite Element Electrochemistry) is a collection of MOOSE-based applications for modeling battery and electrochemical systems. It contains three main applications:

- **eel**: A parallel finite-element code for modeling battery-related physics including mass transport, electrostatics, mechanical deformation, and thermal effects
- **ec_beta**: Electrochemistry beta version (passed testing phase, under development) 
- **ecm_test**: Electro-chemo-mechanics test version (not yet passed testing phase)

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
- `problems/`: Input files (.i) for simulations and test cases
- `test/`: Test framework files
- `doc/`: Documentation and configuration

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