# Phase 2 EDL benchmark with proper AD kernels and SternRobinBC

[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0.0
  xmax = 3.04e-8          # 10 Debye lengths (~30 nm)
  nx = 100
[]

[Variables]
  [./phi]
    initial_condition = 0.0
    scaling = 1.0e-06    # Balanced scaling for electric potential
  [../]
  [./c_plus]
    initial_condition = 10.0   # mol/m^3
    scaling = 1.0e+05    # Balanced scaling for concentrations
  [../]
  [./c_minus]
    initial_condition = 10.0
    scaling = 1.0e+05    # Balanced scaling for concentrations
  [../]
  [./temperature]
    initial_condition = 298.15  # K
    scaling = 1.0
  [../]
[]

[AuxVariables]
  [./charge_density]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./charge_density_calc]
    type = ChargeDensity
    variable = charge_density
    concentrations = 'c_plus c_minus'
    valences = '1.0 -1.0'
    faraday_constant = 96485.33212
  [../]
[]

[Materials]
  [./electrolyte_constants]
    type = EDLProperties
    block = 0
    D_plus = 1.0e-9
    D_minus = 1.0e-9
    permittivity = 6.95104e-10   # epsilon0 * eps_r (78.5)
  [../]
  [./edl_fluxes]
    type = EDLFluxes
    phi = phi
    c_plus = c_plus
    c_minus = c_minus
    permittivity = permittivity
    D_plus = D_plus
    D_minus = D_minus
    faraday_constant = 96485.33212
    gas_constant = 8.314462618
    temperature = 298.15
    z_plus = 1.0
    z_minus = -1.0
  [../]
[]

[Kernels]
  # Poisson equation: -∇·(ε∇φ) = F (z_+ c_+ + z_- c_-)
  [./phi_diffusion]
    type = ADMatDiffusion
    variable = phi
    diffusivity = permittivity
  [../]
  [./phi_charge_source]
    type = ADPoissonChargeRHS
    variable = phi
    c_plus = c_plus
    c_minus = c_minus
    F = 96485.33212
    z_plus = 1.0
    z_minus = -1.0
  [../]

  # Nernst–Planck for c_+ with diffusion and electromigration
  [./cplus_diffusion]
    type = ADMatDiffusion
    variable = c_plus
    diffusivity = D_plus
  [../]
  [./cplus_electromigration]
    type = ADNernstPlanckConvection
    variable = c_plus
    diffusivity = D_plus
    potential = phi
    temperature = temperature
    valence = 1.0
    faraday_constant = 96485.33212
    gas_constant = 8.314462618
  [../]

  # Nernst–Planck for c_- with diffusion and electromigration
  [./cminus_diffusion]
    type = ADMatDiffusion
    variable = c_minus
    diffusivity = D_minus
  [../]
  [./cminus_electromigration]
    type = ADNernstPlanckConvection
    variable = c_minus
    diffusivity = D_minus
    potential = phi
    temperature = temperature
    valence = -1.0
    faraday_constant = 96485.33212
    gas_constant = 8.314462618
  [../]

  # Temperature (constant for this steady-state case)
  [./temperature_dummy]
    type = Reaction
    variable = temperature
    rate = 0.0
  [../]
[]

[BCs]
  # Stern layer Robin BC at electrode (left)
  [./phi_stern]
    type = SternRobinBC
    variable = phi
    boundary = left
    epsilon0 = 8.854187817e-12
    stern_relative_permittivity = 10
    stern_thickness = 2.0e-10
    metal_potential = 1.0e-3
  [../]

  # Reference gauge in bulk
  [./phi_right]
    type = DirichletBC
    variable = phi
    boundary = right
    value = 0.0
  [../]

  # Bulk concentrations at right boundary
  [./cplus_right]
    type = DirichletBC
    variable = c_plus
    boundary = right
    value = 10.0
  [../]
  [./cminus_right]
    type = DirichletBC
    variable = c_minus
    boundary = right
    value = 10.0
  [../]

  # Temperature boundary conditions
  [./temp_left]
    type = DirichletBC
    variable = temperature
    boundary = left
    value = 298.15
  [../]
  [./temp_right]
    type = DirichletBC
    variable = temperature
    boundary = right
    value = 298.15
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  nl_rel_tol = 1e-5    # Relaxed tolerance for demonstration
  nl_abs_tol = 5e-9    # Set absolute tolerance to match balanced residual scale
  nl_max_its = 100     # Increase iteration limit
  l_tol = 1e-10
  l_max_its = 400
  petsc_options_iname = '-pc_type -ksp_type'
  petsc_options_value = 'lu       preonly'
[]

[Debug]
  show_var_residual_norms = true
[]

[Outputs]
  exodus = true
  csv = true
  file_base = rst/comsol_phase2_edl_corrected
[]