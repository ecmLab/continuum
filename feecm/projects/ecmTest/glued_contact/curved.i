[Mesh]
  patch_update_strategy = iteration
  patch_size = 200
  parallel_type = REPLICATED
  [./mesh]
    type = FileMeshGenerator
    file = data/initialDefect2.msh
  [../]
[]

[GlobalParams]
  displacements = 'ux uy'
[]

[Variables]
  [./ux]
  [../]
  [./uy]
  [../]
  [./li_ion_V]
    block = 'blockCeramic blockMetal interLayer'
  [../]
  [./thermal_lm]
    block = 'mech_contact_secondary_subdomain'
  [../]
[]

[AuxVariables]
  [./plastic_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_rate]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./li_metal_conc]
    block = 'interLayer blockMetal'
  [../]
[]

[AuxKernels]
  [./peeq]
    type = ADMaterialRealAux
    variable = plastic_strain
    property = effective_plastic_strain
    execute_on = timestep_end
    block = 'interLayer blockMetal'
  [../]
  [./strain_rate]
    type = ADMaterialRealAux
    variable = strain_rate
    property = plastic_strain_rate
    block = 'interLayer blockMetal'
  [../]
  [./li_metal_conc]
    type = ConstantAux
    value = 0
    variable = li_metal_conc
    block = 'interLayer blockMetal'
  [../]
[]


[Functions]
  [./pressure]
    type = ParsedFunction
    value = 'if (t <= 10, 0.32e-6*t, 3.2e-6)'
  [../]
  [./uy]
    type = ParsedFunction
    value = '0.01*t'
  [../]
  [./vel2]
    type = ParsedFunction
    value = '0.002*exp(0.0001*t)'
  [../]
  [./flux]
    type = ParsedFunction
    value = 'if (t > 40.0, 1e-12, 0.0)'
  [../]
  [./gapk]
    type = PiecewiseLinear
    data_file = data/gap_cond.csv
    format = columns
  [../]
  [./gapk1]
    type = PiecewiseLinear
    x = '-10.0 0 1e-5 1e-4 1 10 100'
    y = '1e-2 1e-2 1e-2 0 0 0 0'
  [../]

[]

[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress'
    use_automatic_differentiation = true
    extra_vector_tags = 'ref'
    block = 'blockMetal blockCeramic interLayer'
    # use_finite_deform_jacobian = true
  [../]
[]
[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'ux uy mech_contact_normal_lm; li_ion_V thermal_lm'
  acceptable_iterations = 2
[]

[Contact]
  [./mech_contact]
    disp_x = ux
    disp_y = uy
    secondary  =  'blockCeramic_top'
    primary =  'blockMetal_bottom'
    normal_smoothing_distance = 0.5
    tangential_tolerance = 0.25
    formulation = mortar
    mesh = mesh
  [../]
[]
[Constraints]
  [./thermal_constraint]
    type = PressureLMConductanceConstraint
    variable = thermal_lm
    secondary_variable = li_ion_V
    secondary_boundary =  'blockCeramic_top'
    secondary_subdomain = 'mech_contact_secondary_subdomain'
    primary_boundary = 'blockMetal_bottom'
    primary_subdomain = 'mech_contact_primary_subdomain'
    k_function = gapk
    k = 1e-8
    pressure = mech_contact_normal_lm
    use_displaced_mesh = true
    compute_lm_residuals = true
  [../]
[]

[Kernels]
  [./li_ion_V]
    type = ADHeatConduction
    block = 'blockCeramic interLayer blockMetal'
    variable = li_ion_V
    use_displaced_mesh = false
  [../]
[]

[Materials]
  [./thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1.0e-8
    block = 'blockCeramic'
  [../]

  [./thermal_conductivity2]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1e-4
    block = 'interLayer blockMetal'
    use_displaced_mesh = true
  [../]


  [./elasticity_tensor_Li]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 7.81e-3
    poissons_ratio = 0.38
    block = 'blockMetal interLayer'
  [../]
  [./elasticity_tensor_Llzo]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 170e-3
    poissons_ratio = 0.3
    block = 'blockCeramic'
  [../]
  [./stress_llzo]
    type = ADComputeFiniteStrainElasticStress
    block = 'blockCeramic'
  [../]

  [./stress_Li]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    block = 'interLayer'
  [../]

  [./plas]
    type = ADIsoTropicHyperViscoSwellingCreep
    hardening_exponent = 2.0
    saturation_resistance = 2.0e-6
    initial_resistance = 0.95e-6
    hardening_modulus = 10.0e-6
    rate_exponent = 0.15
    activation_energy = 37
    gas_constant = 8.314462681e-3 # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = 298
    alpha = '1 0 0'
    intBnd = 'blockCeramic_top'
    concentration = li_metal_conc
    cref = 0
    omega = 2
    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 2.0
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    block = 'interLayer'
  [../]
  [./stress_Li2]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas2'
    block = 'blockMetal'
  [../]

  [./plas2]
    type = ADIsoTropicHyperViscoSwellingCreep
    hardening_exponent = 2.0
    saturation_resistance = 2.0e-6
    initial_resistance = 0.95e-6
    hardening_modulus = 10.0e-6
    rate_exponent = 0.15
    activation_energy = 37
    gas_constant = 8.314462681e-3 # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = 298
    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 2.0
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    block = 'blockMetal'

    alpha = '1 0 0'
    intBnd = 'blockCeramic_top'
    concentration = li_metal_conc
    cref = 0
    omega = 0.0
  [../]
[]

[BCs]
  [./pressure]
    type = PresetVelocity
    variable = uy
    function = vel2
    boundary = 'blockCeramic_bottom'
  [../]

  [./left_right]
    type = ADDirichletBC
    preset =  true
    variable = ux
    boundary = 'blockCeramic_left blockCeramic_bottom blockCeramic_right blockMetal_right blockMetal_left blockMetal_top'
    value = 0
  [../]
  [./bottom]
    type = ADDirichletBC
    preset = true
    variable = uy
    boundary = 'blockMetal_top'
    value = 0
  [../]
  [./li_ion_bottom]
    type = FunctionNeumannBC
    boundary = 'blockCeramic_bottom'
    variable = li_ion_V
    function = flux
    extra_vector_tags = 'ref'
  [../]

  [./li_ion_top]
    type = ADDirichletBC
    boundary = 'blockMetal_top'
    variable = li_ion_V
    value = 0.0
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./matl_ts_min]
    type = MaterialTimeStepPostprocessor
    block = 'interLayer blockMetal'
  [../]
  [./bot_stress]
    type = SideAverageValue
    variable = stress_yy
    boundary = 'blockCeramic_bottom'
  [../]
[]

[Executioner]

  type = Transient
  solve_type = 'NEWTON'
  automatic_scaling = true
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor -snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_force_iteration -pc_factor_shift_type -pc_factor_shift_amount -mat_mffd_err'
  petsc_options_value = 'lu superlu_dist 1     NONZERO               1e-15               1e-5'



  line_search = 'none'


  nl_abs_tol = 3e-11
  nl_rel_tol = 1e-6

  l_max_its = 100
  nl_max_its = 35

  start_time = 0.0
  dt = 0.001
  dtmax = 0.05
  dtmin = 1e-5
  end_time = 40.0
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.001
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 30
    timestep_limiting_postprocessor = matl_ts_min
  [../]

[] # Executioner

[Outputs]
  [./out]
    type = Exodus
    sync_only = false
    file_base = rst/curved
    interval = 3
  [../]
  [./csv]
    type = CSV
    file_base = rst/curved
  [../]
  [./check]
    type = Checkpoint
    file_base = check/curved
    num_files = 5
    start_time = 0.0
    sync_only = false
    sync_times = '15.0 17.0 20.0 25.0 30.0 32.0 35.0 40.0'
    interval = 50
  [../]
[] # Outputs
