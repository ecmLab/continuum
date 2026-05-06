# Example demonstrating the NonDimensionalParameters material class
# Copyright 2025, CEWLAB, All Rights Reserved

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  xmax = 1.0
  ymax = 1.0
[]

[Variables]
  [concentration]
    order = FIRST
    family = LAGRANGE
  []
  [potential]
    order = FIRST
    family = LAGRANGE
  []
[]

[Materials]
  # NonDimensionalParameters material with Li-ion battery characteristic scales
  [nondim_params]
    type = NonDimensionalParameters
    # Li-ion battery characteristic scales
    c0 = 29400.0    # Li concentration in Li metal [mol/m³]
    T0 = 298.0      # Reference temperature [K]
    D0 = 1.0e-14    # Li diffusivity in graphite [m²/s]
    j0 = 1.0        # Exchange current density [A/m²]
    Omega0 = 1.304e-5  # Li molar volume [m³/mol]
    F = 96485.0     # Faraday constant [C/mol]
    R = 8.314       # Gas constant [J/mol/K]
  []
  
  # Example material showing how to access the characteristic scales
  [example_diffusivity]
    type = GenericConstantMaterial
    prop_names = 'diffusivity'
    prop_values = '1.0'  # Will be made dimensionless using D0
  []
[]

[Kernels]
  # Simple diffusion kernel for demonstration
  [diffusion_c]
    type = Diffusion
    variable = concentration
  []
  [diffusion_phi]
    type = Diffusion
    variable = potential
  []
[]

[BCs]
  [left_c]
    type = DirichletBC
    variable = concentration
    boundary = left
    value = 1.0
  []
  [right_c]
    type = DirichletBC
    variable = concentration
    boundary = right
    value = 0.0
  []
  [left_phi]
    type = DirichletBC
    variable = potential
    boundary = left
    value = 1.0
  []
  [right_phi]
    type = DirichletBC
    variable = potential
    boundary = right
    value = 0.0
  []
[]

[Postprocessors]
  # Access characteristic scales through material properties
  [L0_value]
    type = ElementAverageMaterialProperty
    mat_prop = L0
  []
  [t0_value]
    type = ElementAverageMaterialProperty
    mat_prop = t0
  []
  [phi0_value]
    type = ElementAverageMaterialProperty
    mat_prop = phi0
  []
  [sigma0_value]
    type = ElementAverageMaterialProperty
    mat_prop = sigma0
  []
  [kappa_value]
    type = ElementAverageMaterialProperty
    mat_prop = kappa
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-8
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
[]