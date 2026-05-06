# Test for full cell
# Geometry consists of
# 1) Anode layer with given OCV curve = nmc_equilibrium_potential.csv
# 2) Solid Electrolyte
# 3) Anode at 0 V
# No mechanics included here for simplicity
# The goal is to just check to see of the whole electro-chemistry works
# All the GapDisplacementConductanceConstraint is just Linear Butler Volmer Kinetics
# All voltages are in mV => 3V = 3000 mV
# --- Issues ----
# 1) Changing the charge transfer resistance changes the voltage profile
# 2) Convergence changes when changing R_ct
# Note k in GapDisplacementConductanceConstraint is 1/Rct
[Mesh]
  [./mesh]
    type = FileMeshGenerator
    file = 'test_anode.msh'
  [../]
  [./primary_block]
    type = LowerDBlockFromSidesetGenerator
    input = mesh
    sidesets = 'Electrolyte_bottom'
    new_block_name = 'primary_contact_block'
  [../]
  [./secondary_block]
    type = LowerDBlockFromSidesetGenerator
    input = primary_block
    sidesets = 'Anode_top'
    new_block_name = 'secondary_contact_block'
  [../]

[]
[GlobalParams]
  displacements = 'ux uy'
[]

[Variables]
  [./V]
    block = 'Anode Electrolyte'
  [../]
  [./thermal_lm]
    block = 'secondary_contact_block'
  [../]
  [./li_metal_conc]
    initial_condition = 0.0
    block = 'Anode'
  [../]
  [./ux]
    block = 'Anode Electrolyte'
  [../]
  [./uy]
    block = 'Anode Electrolyte'
  [../]
  [./normal_lm_x_anode]
    block = 'secondary_contact_block'
  [../]

  [./normal_lm_y_anode]
    block = 'secondary_contact_block'
  [../]

[]

[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress hydrostatic_stress'
    use_automatic_differentiation = true
    block = 'Anode Electrolyte'
    # extra_vector_tags = 'ref'
    # use_finite_deform_jacobian = true
  [../]
[]


[Constraints]
  [./thermal_constraint]
    type = GapDisplacementConductanceConstraint
    variable = thermal_lm
    secondary_variable = V
    secondary_boundary =  'Anode_top'
    secondary_subdomain = 'secondary_contact_block'
    primary_boundary = 'Electrolyte_bottom'
    primary_subdomain = 'primary_contact_block'
    # k_function = gapk
    k = 1e-4
    # use_displaced_mesh = true
    displacements = 'ux uy'
    compute_lm_residuals = true
    include_equilibrium_potential = false
    # extra_vector_tags = 'ref'
  [../]

  [./anode_conc]
    type = ScaledBCConstraint
    variable = thermal_lm
    secondary_variable = li_metal_conc
    primary_variable = V
    primary = false
    secondary_boundary =  'Anode_bottom'
    secondary_subdomain = 'secondary_contact_block'
    primary_boundary = 'Electrolyte_top'
    primary_subdomain = 'primary_contact_block'
    scale = 1.036426e-2
    extra_vector_tags = 'ref'
  [../]

  [./anode_mech_x]
    type = EqualValueConstraint
    variable = normal_lm_x_anode
    primary_variable = ux
    secondary_variable = ux
    secondary_boundary =  'Anode_bottom'
    secondary_subdomain = 'secondary_contact_block'
    primary_boundary = 'Electrolyte_top'
    primary_subdomain = 'primary_contact_block'
    use_displaced_mesh = true
  [../]
  [./anode_mech_y]
    type = EqualValueConstraint
    variable = normal_lm_y_anode
    primary_variable = uy
    secondary_variable = uy
    secondary_boundary =  'Anode_bottom'
    secondary_subdomain = 'secondary_contact_block'
    primary_boundary = 'Electrolyte_top'
    primary_subdomain = 'primary_contact_block'
    use_displaced_mesh = true
  [../]


[]


[AuxVariables]
  # [./li_metal_conc]
  #   block = 'Anode'
  # [../]


  [./flux_x]
    order = FIRST
    family = MONOMIAL
    block = 'Anode Electrolyte'
  [../]
  [./flux_y]
    order = FIRST
    family = MONOMIAL
    block = 'Anode Electrolyte'
  [../]
  [./flux_z]
    order = FIRST
    family = MONOMIAL
    block = 'Anode Electrolyte'
  [../]
  [./bndliflux]
    order = CONSTANT
    family = MONOMIAL
    block = 'Anode Electrolyte'
  [../]
  # [./Eq_pot]
  #   order = CONSTANT
  #   family = MONOMIAL
  #   block = 'Anode'
  # [../]
[]

[AuxKernels]
  # [./conc]
  #   type = FunctionAux
  #   function = conc
  #   variable = li_metal_conc
  #   block = 'Anode'
  # [../]
  # [./ux]
  #   type = ConstantAux
  #   value = 0
  #   block = 'Anode Electrolyte'
  #   variable = ux
  # [../]
  # [./uy]
  #   type = ConstantAux
  #   value = 0
  #   block = 'Anode Electrolyte'
  #   variable = uy
  # [../]

  # [./Eq_pot]
  #   type = ADMaterialRealAux
  #   variable = Eq_pot
  #   property = equilibrium_potential
  #   block = 'Anode'
  # [../]

  [./li_ion_flux_x]
    type = ADDiffusionFluxAux
    variable = flux_x
    component = x
    diffusion_variable = V
    diffusivity = thermal_conductivity
    block = 'Anode Electrolyte'
  [../]

  [./li_ion_flux_y]
    type = ADDiffusionFluxAux
    variable = flux_y
    component = y
    diffusion_variable = V
    diffusivity = thermal_conductivity
    block = 'Anode Electrolyte'
  [../]
  [./bnd_li_flux]
    type = DiffusionFluxNormalToBoundaryAux
    variable = bndliflux
    boundary = 'Electrolyte_bottom Anode_top'
    diffusion_variable = V
    diffusivity = thermal_conductivity
  [../]
[]


[Kernels]
  [./li_ion_V]
    type = ADHeatConduction
    variable = V
    use_displaced_mesh = false
    block = 'Anode Electrolyte'
  [../]
  [./li_metal_anode]
    type = ADChemoMechanoAnsioDiffusion
    block = 'Anode'
    variable = li_metal_conc
    diffusivity = diffusivity
    use_displaced_mesh = false
  [../]

  [./li_metal_dt]
    type = ADTimeDerivative
    block = 'Anode'
    variable = li_metal_conc
    use_displaced_mesh = false
  [../]


[]

[Materials]

  [./thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 100
    block = 'Anode'
  [../]
  [./thermal_conductivity2]
    type = ADHeatConductionMaterial
    thermal_conductivity = 0.1
    block = 'Electrolyte'
  [../]
  [./diffusivity_Li_anode]
    type = ADDiffusionAlongPrincipalDirectionsMaterial
    diffusivity_vector = '1e8 0 0'
    block = 'Anode'
  [../]

  [./elasticity_tensor_Li]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 7.8e3
    poissons_ratio = 0.3
    block = 'Anode'
  [../]
  [./elasticity_tensor_Llzo]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 170e3
    poissons_ratio = 0.25
    block = 'Electrolyte'
  [../]
  [./stress_llzo]
    type = ADComputeFiniteStrainElasticStress
    block = 'Electrolyte'
  [../]
  [./stress_Li]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    # perform_finite_strain_rotations = true
    block = 'Anode'
  [../]
  [./plas]
    type = ADIsoTropicHyperViscoSwellingCreep
    absolute_tolerance = 1e-8
    relative_tolerance = 1e-6
    hardening_exponent = 2.0
    saturation_resistance = 2.0
    initial_resistance = 0.95
    hardening_modulus = 10.0
    rate_exponent = 0.15
    # reference_strain_rate = 0.05
    activation_energy = 37e3
    gas_constant = 8.314462681 # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = 298

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 2.0
    omega = 13
    alpha = '1 0 0'
    concentration = li_metal_conc
    cref = 0
    block = 'Anode'
    intBnd = 'Electrolyte_top'
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    enable = true
  [../]


[]

[BCs]

  [./current]
    type = ADFunctionNeumannBC
    boundary = 'Electrolyte_current'
    # function = 'if (x > 4000, 0, 20e-3)'
    # value = 20e-3
    function = flux
    variable = V
    extra_vector_tags = 'ref'
  [../]
  [./OV]
    type = ADDirichletBC
    boundary = 'Anode_bottom'
    value = 0.0
    variable = V
  [../]
  [./ux]
    type = ADDirichletBC
    boundary = 'Electrolyte_top Electrolyte_left Electrolyte_right'
    variable = ux
    value = 0
  [../]
  [./uy]
    type = ADDirichletBC
    boundary = 'Electrolyte_top Electrolyte_left Electrolyte_right'
    variable = uy
    value = 0
  [../]
  [./pre]
    type = ADPressure
    boundary = 'Anode_bottom'
    variable = uy
    component = 1
    function = pressure
  [../]
[]

[Functions]
  [./flux]
    type = ParsedFunction
    value  = 'if (t <= 10, 2e-3*t, 20e-3)'
  [../]
  [./pressure]
    type = ParsedFunction
    value = '1e-1'
  [../]
[]
[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./Voltage_Anode]
    type = SideAverageValue
    boundary = 'Electrolyte_top'
    variable = V
  [../]

  # [./eq_pot]
  #   type = SideAverageValue
  #   boundary = 'Anode_bottom'
  #   variable = Eq_pot
  # [../]
  # [./bot_conc]
  #   type = SideAverageValue
  #   boundary = 'Anode_bottom'
  #   variable = li_metal_conc
  # [../]
[]
# [UserObjects]
#   [./vcutoff]
#     type = Terminator
#     expression = 'Voltage_Anode < 2750'
#     execute_on = 'TIMESTEP_END'
#   [../]
#   [./max_conc]
#     type = Terminator
#     expression = 'bot_conc > 4.77e-2'
#     execute_on = 'TIMESTEP_END'
#   [../]
# []

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  automatic_scaling = true
  compute_scaling_once = true
  # petsc_options_iname = '-pc_type -pc_mat_solver_package'
  # petsc_options_value = 'lu superlu_dist'
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1     NONZERO               1e-15               '
  dt = 200
  # line_search = basic
  # num_steps = 20
  nl_max_its = 25
  nl_abs_tol = 1e-13
  nl_rel_tol = 1e-3
  dtmax = 200
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 15
    # timestep_limiting_postprocessor = matl_ts_min
  [../]
  # end_time = 10000

  num_steps = 2
  snesmf_reuse_base = true
  # scaling_group_variables = 'V thermal_lm'

[]
[Problem]
  type = ReferenceResidualProblem
  # solution_variables = 'ux uy normal_lm thermal_lm li_ion_V thermal_lm2 li_metal_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'V thermal_lm; ux uy normal_lm_x_anode normal_lm_y_anode'
  acceptable_iterations = 2
  # coord_type = RZ
[]
[Outputs]
  exodus = true
  csv = true
[]
