# Test for Si to compare to Bower and Guduru (2014)
# Material Properties
# Molar density of Si = c_0 = 7.874e4 mol/m^3 = 7.874e-2
# Modulus of Si = 100 GPa = 100e3
# poissons_ratio of Si = 0.26
# Vol. expansion of Si (normalized conc ) = 0.7
#    -> Vol. expansion of Si -> 0.7 * 7.874e-2 = omega
# epsilon_0 (strain rate Si) = 0.6e-9 /s
# initial yield stress of Si = 0.12 GPa = 0.12e3
# Stress exponent for plastic flow in Si = 4 -> 0.25
# initial conc of Li in Si = 0.0078 = 0.0078 * 7.87e-2
#applied current density = 0.012 mA/cm^2 = 1.2e-5
# exchange_current_density = 0.001 mA/cm^2 = 1e-6
[Mesh]
  [./Si]
    type = GeneratedMeshGenerator
    nx = 2
    ny = 25
    xmin = 0
    xmax = 0.01
    ymin = 0
    ymax = 0.2
    elem_type = QUAD4
    dim = 2
  [../]
[]
# [GlobalParams]
#   displacements = 'ux uy'
# []

[Variables]
  [./V]
  [../]
  [./li_metal_conc]
  [../]
[]

[Functions]
  [./reaction_rate]
    type = PiecewiseLinear
    x = '0.0 4.857'
    y = '780.0 0.0'
  [../]
[]

[AuxVariables]
  [./V0]
  [../]
  [./flux_x]
    order = FIRST
    family = MONOMIAL
  [../]
  [./flux_y]
    order = FIRST
    family = MONOMIAL
  [../]
  [./flux_z]
    order = FIRST
    family = MONOMIAL
  [../]
  [./bndliflux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Eq_pot]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./bnd_li_flux]
    type = DiffusionFluxNormalToBoundaryAux
    variable = bndliflux
    boundary = 'top'
    diffusion_variable = V
    diffusivity = thermal_conductivity
  [../]
  [./Eq_pot]
    type = ADMaterialRealAux
    variable = Eq_pot
    property = equilibrium_potential
  [../]

  [./V0]
    type = ConstantAux
    value = 780.0
    variable = V0
  [../]
  [./li_ion_flux_x]
    type = ADDiffusionFluxAux
    variable = flux_x
    component = x
    diffusion_variable = V
    diffusivity = thermal_conductivity
  [../]

  [./li_ion_flux_y]
    type = ADDiffusionFluxAux
    variable = flux_y
    component = y
    diffusion_variable = V
    diffusivity = thermal_conductivity
  [../]
  [./li_ion_flux_z]
    type = ADDiffusionFluxAux
    variable = flux_z
    component = z
    diffusion_variable = V
    diffusivity = thermal_conductivity
  [../]
[]


[Kernels]
  [./li_ion_V]
    type = ADHeatConduction
    variable = V
    use_displaced_mesh = false
  [../]
  [./li_metal2]
    type = ADMatDiffusion
    variable = li_metal_conc
    diffusivity = diffusivity
    use_displaced_mesh = false
  [../]
  [./li_metal_dt]
    type = ADTimeDerivative
    variable = li_metal_conc
    use_displaced_mesh = false
  [../]

[]

[Materials]
  [./thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1
  [../]
  [./diffusivity_Li]
    type = ADGenericConstantMaterial
    prop_names = 'diffusivity'
    prop_values = '3e-7'
  [../]
  [./equilibrium_potential]
    type = ADComputeEquilibriumPotential
    R = 8.31446
    faraday = 96.4853329
    temperature = 298
    cref = 7.874e-2
    concentration = li_metal_conc
    include_conc = false
    include_reaction_rate = true
    # reaction_rate = 780.0
    reaction_rate_function = reaction_rate
    include_mechanical_effects = false
  [../]

[]

[BCs]
  [./current]
    type = ADButlerVolmerBC
    current_density = 1.2e-5
    exchange_current_density = 1e-6
    faraday = 96.4853
    Temperature = 298
    R = 8.314462681
    variable = V
    boundary = top
    # extra_vector_tags = 'ref'
  [../]
  [./conc]
    type = ADNeumannBC
    variable = li_metal_conc
    # v = bndliflux
    # scale = -1.036428e-2
    value = 1.24371375e-7
    boundary = top
    # extra_vector_tags = 'ref'
  [../]

[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  automatic_scaling = true
  compute_scaling_once = false
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration'# -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1'#     NONZERO               1e-20               '
  dt = 200
  # num_steps = 20
  nl_max_its = 35
  nl_abs_tol = 5e-11
  nl_rel_tol = 1e-6
  dtmax = 50
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 50
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 40
    # timestep_limiting_postprocessor = matl_ts_min
  [../]
  end_time = 8000
  snesmf_reuse_base = false
  # scaling_group_variables = 'ux uy; V li_metal_conc'

[]

[Outputs]
  exodus = true
[]
