# elem = QUAD4
# order = FIRST
[Mesh]
  patch_size = 80
  patch_update_strategy = auto
  parallel_type = REPLICATED
  [./mesh]
    type = FileMeshGenerator
    file = data/cosine_symmetric2.msh
  [../]
[]

[GlobalParams]
  displacements = 'ux uy'
[]

[Variables]
  [./ux]
    block = 'blockCeramic blockMetal interLayer'
  [../]
  [./uy]
    block = 'blockCeramic blockMetal interLayer'
  [../]
  []

[AuxVariables]
  # [./plastic_strain]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./strain_rate]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  [./li_metal_conc]
  [../]
[]
[AuxKernels]
  [./conc]
    type = ConstantAux
    value = 0.0
    variable = li_metal_conc
  [../]
  # [./peeq]
  #   type = ADMaterialRealAux
  #   variable = plastic_strain
  #   property = effective_plastic_strain
  #   execute_on = timestep_end
  #   block = 'interLayer blockMetal'
  # [../]
  # [./strain_rate]
  #   type = ADMaterialRealAux
  #   variable = strain_rate
  #   property = plastic_strain_rate
  #   block = 'interLayer blockMetal'
  # [../]
[]
[Functions]
  [./pressure]
    type = PiecewiseLinear
    data_file = 'pressure_time1.csv'
    format = columns
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress hydrostatic_stress'
    use_automatic_differentiation = true
    block = 'blockCeramic interLayer blockMetal'
    extra_vector_tags = 'ref'
    # use_finite_deform_jacobian = true
  [../]
[]

[Problem]
  type = ReferenceResidualProblem
  # solution_variables = 'ux uy normal_lm thermal_lm li_ion_V thermal_lm2 li_metal_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'ux uy mech_normal_lm'
  acceptable_iterations = 2
[]

[Contact]
  [./mech]
    primary = 'blockCeramic_top'
    secondary = 'blockMetal_bottom'
    model = frictionless
    formulation = mortar
    mesh = mesh
    # use_dual = true
    c_normal = 1e-3
    normal_smoothing_distance = 0.2
    tangential_tolerance = 0.25
  [../]

[]

[Materials]
  [./elasticity_tensor_Li]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 7.8e-3
    poissons_ratio = 0.3
    block = 'interLayer blockMetal'
  [../]
  [./elasticity_tensor_Llzo]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 170e-3
    poissons_ratio = 0.25
    block = 'blockCeramic'
  [../]
  [./stress_llzo]
    type = ADComputeFiniteStrainElasticStress
    block = 'blockCeramic'
  [../]

  [./stress_Li]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    # perform_finite_strain_rotations = true
    block = 'interLayer'
  [../]

  [./stress_Li2]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas2'
    # perform_finite_strain_rotations = true
    block = 'blockMetal'
  [../]
  [./plas]
    type = ADIsoTropicHyperViscoCreep
    # absolute_tolerance = 1e-6
    # relative_tolerance = 1e-06
    hardening_exponent = 2.0
    saturation_resistance = 2.0e-6
    initial_resistance = 0.95e-6
    hardening_modulus = 10.0e-6
    rate_exponent = 0.15
    # reference_strain_rate = 0.05
    activation_energy = 37
    gas_constant = 8.314462681e-3 # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = 333

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 0.1
    omega = 2
    # alpha = '1 0 0'
    # concentration = 0.0
    # cref = 0
    block = 'interLayer'
    # intBnd = 'blockMetal_bottom'
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    enable = true
  [../]

  [./plas2]
    type = ADIsoTropicHyperViscoCreep
    hardening_exponent = 2.0
    saturation_resistance = 2.0e-6
    initial_resistance = 0.95e-6
    hardening_modulus = 10.0e-6
    rate_exponent = 0.15
    activation_energy = 37
    gas_constant = 8.314462681e-3 # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = 333

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 0.1
    # omega = 2e-6
    # alpha = '1 0 0'
    # intBnd = 'blockMetal_top'
    # concentration = 0.0
    # cref = 0
    block = 'blockMetal'
  [../]
[]

[BCs]
  [./Li_top_y]
    type = ADDirichletBC
    variable = uy
    boundary = 'blockMetal_top'
    value = 0.0
  [../]

  [./Li_top_x]
    type = ADDirichletBC
    boundary = 'blockMetal_top'
    variable = ux
    value = 0.0
  [../]

  [./left_right]
    type = ADDirichletBC
    variable = ux
    boundary = 'blockCeramic_left blockCeramic_right blockMetal_left blockMetal_right'
    value = 0
  [../]

  [./top_pressure]
    type = PresetVelocity
    boundary = 'blockCeramic_bottom'
    variable = uy
    # component = 1
    velocity = 1e-6
    # function = pressure
    function = 'if (t < 5, 1, 0)'
    # use_displaced_mesh = true
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

  #Preconditioned JFNK (default)
  solve_type = 'NEWTON'
  automatic_scaling = true
  compute_scaling_once = false
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor'
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1     NONZERO               1e-20               '
  line_search = contact
  l_max_its = 50
  nl_max_its = 25
  nl_abs_tol = 3e-15
  nl_rel_tol = 1e-10
  start_time = 0.0
  dt = 1.0
  dtmax = 2.0
  dtmin = 1e-5
  end_time = 500.0
  # num_steps = 10
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.025
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 40
    # timestep_limiting_postprocessor = matl_ts_min
  [../]

[] # Executioner

[Outputs]
  [./out]
    type = Exodus
  [../]
  [./csv]
    type = CSV
  [../]
[] # Outputs

[Debug]
  show_var_residual_norms = true
  show_material_props = true
[]
