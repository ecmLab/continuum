# Example demonstrating non-dimensional diffusion using the TECM framework
# Copyright 2025, CEWLAB, All Rights Reserved

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
  xmax = 1.0  # This will be interpreted as x̃_max = 1 (dimensionless)
  ymax = 1.0  # This will be interpreted as ỹ_max = 1 (dimensionless)
[]

[Variables]
  # Dimensionless concentration c̃ = c/c₀
  [concentration_nondim]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.1
  []
[]

[Materials]
  # Characteristic scales for Li-ion battery system (evaluated first)
  [nondim_params]
    type = NonDimensionalParameters
    # Li metal + Li6PS5Cl + NMC battery system characteristic scales
    c0 = 41528.0       # Li concentration in Li6PS5Cl solid electrolyte [mol/m³] (primary scale)
    T0 = 298.0         # Reference temperature [K]
    D0 = 2.5e-12       # Li diffusivity in Li6PS5Cl solid electrolyte [m²/s] (reference)
    j0 = 1.0           # Exchange current density [A/m²]
    Omega0 = 1.304e-5  # Li molar volume [m³/mol]
    F = 96485.0        # Faraday constant [C/mol]
    R = 8.314          # Gas constant [J/mol/K]
    # This gives: σ₀ = F²c₀D₀/(RT₀) ≈ 0.39 S/m (Nernst-Einstein)
    #            L₀ = Fc₀D₀/j₀ ≈ 10.0 mm
    #            t₀ = L₀²/D₀ ≈ 1.27 years
    outputs = exodus
    output_properties = 'L0 t0 phi0 mu0 sigma0'
  []
  
  # Different diffusivities to demonstrate scaling (evaluated second)
  [diffusivity_electrolyte]
    type = NonDimensionalDiffusivity
    D_dimensional = 1.0e-13  # Li diffusivity in liquid electrolyte [m²/s]
    # This gives D̃ = 1.0e-13 / 2.5e-12 = 0.04
  []
[]

[Kernels]
  # Non-dimensional diffusion: -∇̃ · (D̃ ∇̃c̃) = 0
  [diffusion_nondim]
    type = NonDimensionalDiffusion
    variable = concentration_nondim
    diffusivity_nondim = diffusivity_nondim
  []
[]

[BCs]
  # Boundary conditions in dimensionless form
  [left_high]
    type = DirichletBC
    variable = concentration_nondim
    boundary = left
    value = 1.0  # c̃ = 1 means c = c₀
  []
  [right_low]
    type = DirichletBC
    variable = concentration_nondim
    boundary = right
    value = 0.1  # c̃ = 0.1 means c = 0.1c₀
  []
  [top_insulated]
    type = NeumannBC
    variable = concentration_nondim
    boundary = top
    value = 0.0  # No flux
  []
  [bottom_insulated]
    type = NeumannBC
    variable = concentration_nondim
    boundary = bottom
    value = 0.0  # No flux
  []
[]

[Postprocessors]
  # Monitor characteristic scales
  [L0_microns]
    type = ElementAverageMaterialProperty
    mat_prop = L0
  []
  [t0_hours]
    type = ElementAverageMaterialProperty
    mat_prop = t0
  []
  [D_nondim_avg]
    type = ADElementAverageMaterialProperty
    mat_prop = diffusivity_nondim
  []
  
  # Monitor solution behavior
  [c_nondim_max]
    type = NodalExtremeValue
    variable = concentration_nondim
    value_type = max
  []
  [c_nondim_min]
    type = NodalExtremeValue
    variable = concentration_nondim
    value_type = min
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-10
  l_max_its = 100
  nl_max_its = 20
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
[]