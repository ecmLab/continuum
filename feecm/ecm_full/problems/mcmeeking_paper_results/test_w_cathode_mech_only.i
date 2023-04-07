# Test to check results of mcmeeking paper
# Test overpotential on the surface etc.

current_density = 10e-3
k_anode = 1e-3
faraday = 96.4853329
temperature = 298
gas_constant = 8.31446
c_ref = 0
c_ref_pore = 2.5925e-2
k_pore = 1e-3
c_init_Li = 0
c_init_C = 2.5925e-2
c_max_C = 3.05e-2
pressure = 1e-4
end_time = 150
[Mesh]
  construct_node_list_from_side_list = true
  patch_update_strategy=iteration
  [./mesh]
    type = FileMeshGenerator
    file = 'data/test_w_cathode2.e'
  [../]
  # ------------ Interfaces ---------------
  # --- interface between interlayer and ceramic
  [./secondary_interLayer_block]
    type = LowerDBlockFromSidesetGenerator
    input = mesh
    sidesets = 'interLayer_bottom'
    new_block_name = 'interLayer_s_block'
  [../]
  [./primary_interLayer_block]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_interLayer_block
    sidesets = 'blockCeramic_top'
    new_block_name = 'interLayer_p_block'
  [../]
  # --- interface between pore and ceramic
  [./secondary_Pore_block]
    type = LowerDBlockFromSidesetGenerator
    input = primary_interLayer_block
    sidesets = 'blockPore_bottom'
    new_block_name = 'blockPore_s_block'
  [../]
  [./primary_Pore_block]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_Pore_block
    sidesets = 'blockCeramic_top'
    new_block_name = 'blockPore_p_block'
  [../]
  # --- interface between ceramic and cathode
  [./secondary_Cathode_block]
    type = LowerDBlockFromSidesetGenerator
    input = primary_Pore_block
    sidesets = 'blockCeramic_bottom'
    new_block_name = 'blockCathode_s_block'
  [../]
  [./primary_Cathode_block]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_Cathode_block
    sidesets = 'blockCathode_top'
    new_block_name = 'blockCathode_p_block'
  [../]

[]

[GlobalParams]
  displacements = 'ux uy'
[]
[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress hydrostatic_stress'
    use_automatic_differentiation = true
    block = 'blockCeramic interLayer blockPore blockMetal blockCathode '
    extra_vector_tags = 'ref'
  [../]
[]

[Contact]
  [./metal_interlayer]
    formulation = mortar
    model = frictionless
    primary = 'interLayer_top'
    secondary ='blockMetal_bottom'
    normal_smoothing_distance = 0.5
    tangential_tolerance = 0.25
    mesh = primary_Cathode_block
  [../]
  [./metal_pore]
    formulation = mortar
    model = frictionless
    secondary =  'blockMetal_right'
    primary = 'blockPore_left'
    normal_smoothing_distance = 0.5
    tangential_tolerance = 0.25
    mesh = metal_interlayer_secondary_subdomain_generator
  [../]
  [./interlayer_pore]
    formulation = mortar
    model = frictionless
    secondary =  'interLayer_right'
    primary  ='blockPore_left'
    normal_smoothing_distance = 0.5
    tangential_tolerance = 0.25
    mesh = metal_pore_secondary_subdomain_generator
  [../]
[]

[Variables]
  [./ux]
    block = 'blockCeramic interLayer blockPore blockMetal blockCathode '
  [../]
  [./uy]
    block = 'blockCeramic interLayer blockPore blockMetal blockCathode '
  [../]

  [./normal_lm_x_interlayer_ceramic]
    # Variables for interlayer_ceramic_block
    block = 'interLayer_s_block'
  [../]
  [./normal_lm_y_interlayer_ceramic]
    block = 'interLayer_s_block'
  [../]
  # Variables for interlayer_ceramic_block
  [./normal_lm_x_pore_ceramic]
    block = 'blockPore_s_block'
  [../]
  [./normal_lm_y_pore_ceramic]
    block = 'blockPore_s_block'
  [../]
  # Variables for interlayer_ceramic_block
  [./normal_lm_x_cathode_ceramic]
    block = 'blockCathode_s_block'
  [../]
  [./normal_lm_y_cathode_ceramic]
    block = 'blockCathode_s_block'
  [../]
[]

[Constraints]
  # ---- Cathode Ceramic
  [./Cathode_ceramic_x]
    type = EqualValueConstraint
    primary_variable = ux
    secondary_variable = ux
    secondary_boundary = 'blockCeramic_bottom'
    secondary_subdomain = 'blockCathode_s_block'
    primary_boundary = 'blockCathode_top'
    primary_subdomain = 'blockCathode_p_block'
    variable = normal_lm_x_cathode_ceramic
    extra_vector_tags = 'ref'
    use_displaced_mesh = true
  [../]
  [./Cathode_ceramic_y]
    type = EqualValueConstraint
    primary_variable = uy
    secondary_variable = uy
    secondary_boundary = 'blockCeramic_bottom'
    secondary_subdomain = 'blockCathode_s_block'
    primary_boundary = 'blockCathode_top'
    primary_subdomain = 'blockCathode_p_block'
    variable = normal_lm_y_cathode_ceramic
    extra_vector_tags = 'ref'
    use_displaced_mesh = true
  [../]
  # ---- Pore Ceramic
  [./Pore_ceramic_x]
    type = EqualValueConstraint
    primary_variable = ux
    secondary_variable = ux
    secondary_boundary = 'blockPore_bottom'
    secondary_subdomain = 'blockPore_s_block'
    primary_boundary = 'blockCeramic_top'
    primary_subdomain = 'blockPore_p_block'
    variable = normal_lm_x_pore_ceramic
    extra_vector_tags = 'ref'
    use_displaced_mesh = true
  [../]
  [./Pore_ceramic_y]
    type = EqualValueConstraint
    primary_variable = uy
    secondary_variable = uy
    secondary_boundary = 'blockPore_bottom'
    secondary_subdomain = 'blockPore_s_block'
    primary_boundary = 'blockCeramic_top'
    primary_subdomain = 'blockPore_p_block'
    variable = normal_lm_y_pore_ceramic
    extra_vector_tags = 'ref'
    use_displaced_mesh = true
  [../]

  # ---- Interlayer Ceramic
  [./interLayer_ceramic_x]
    type = EqualValueConstraint
    primary_variable = ux
    secondary_variable = ux
    secondary_boundary = 'interLayer_bottom'
    secondary_subdomain = 'interLayer_s_block'
    primary_boundary = 'blockCeramic_top'
    primary_subdomain = 'interLayer_p_block'
    variable = normal_lm_x_interlayer_ceramic
    extra_vector_tags = 'ref'
    use_displaced_mesh = true
  [../]
  [./interLayer_ceramic_y]
    type = EqualValueConstraint
    primary_variable = uy
    secondary_variable = uy
    secondary_boundary = 'interLayer_bottom'
    secondary_subdomain = 'interLayer_s_block'
    primary_boundary = 'blockCeramic_top'
    primary_subdomain = 'interLayer_p_block'
    variable = normal_lm_y_interlayer_ceramic
    extra_vector_tags = 'ref'
    use_displaced_mesh = true
  [../]
[]

### --- for mech only problem give aux var for conc
[AuxVariables]
  [./li_metal_conc]
    block = 'interLayer'
  [../]
[]
[AuxKernels]
  [./li_metal_conc]
    type = ConstantAux
    value = 0
    variable = li_metal_conc
  [../]
[]

[Materials]
  # --- Mechanical Elastic Properties
  # --- Cathode mechanical properties
  [./elasticity_Cathode]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 100e3
    poissons_ratio = 0.25
    block = 'blockCathode'
  [../]
  # --- Pore mechanical properties
  [./elasticity_Pore]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 130e3
    poissons_ratio = 0.25
    block = 'blockPore'
  [../]
  # --- Ceramic mechanical properties
  [./elasticity_Ceramic]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 170e3
    poissons_ratio = 0.25
    block = 'blockCeramic'
  [../]
  # --- Li metal mechanical properties
  [./elasticity_Li]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 7.8e3
    poissons_ratio = 0.31
    block = 'blockMetal interLayer'
  [../]
  # ---- Stress calculators
  [./stress_elastic]
    type = ADComputeFiniteStrainElasticStress
    block = 'blockCathode blockCeramic blockPore'
  [../]
  # --- Interlayer Stress calculator
  [./stress_interlayer]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas_interlayer'
    block = 'interLayer'
  [../]
  # --- Interlayer Stress calculator
  [./stress_Li_metal]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas_Li'
    block = 'blockMetal'
  [../]
  # --- interLayer Plasticity
  [./plas_interlayer]
    type = ADIsoTropicHyperViscoSwellingCreep
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
    alpha = '1 0 0'
    concentration = li_metal_conc
    cref = 0
    block = 'interLayer'
    intBnd = 'blockMetal_bottom'
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    enable = true
  [../]
  # --- Li metal plasticity
  [./plas_Li]
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
    block = 'blockMetal'
  [../]
[]

[Functions]
  [./pressure]
    type = ParsedFunction
    value = '1e-2*t'
  [../]
[]

[BCs]
  [./bottom_y]
    type = ADDirichletBC
    variable = uy
    boundary = 'blockCathode_bottom'
    preset = true
    value = 0
  [../]
  [./left_symmetry]
    type = ADDirichletBC
    variable = ux
    boundary = 'blockCathode_left blockCeramic_left interLayer_left blockMetal_left blockPore_right'
    value = 0
  [../]
  [./top_pressure]
    type = ADPressure
    boundary = 'blockMetal_top'
    variable = uy
    component = 1
    function = pressure
    use_displaced_mesh = true
    extra_vector_tags = 'ref'
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
  line_search = none
  l_max_its = 50
  nl_max_its = 25
  nl_abs_tol = 1e-16
  nl_rel_tol = 1e-4
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

[] # Executione

[Problem]
  type = ReferenceResidualProblem
  # solution_variables = 'ux uy normal_lm thermal_lm li_ion_V thermal_lm2 li_metal_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = "ux uy normal_lm_x_interlayer_ceramic normal_lm_y_interlayer_ceramic normal_lm_x_pore_ceramic normal_lm_y_pore_ceramic normal_lm_x_cathode_ceramic
                            normal_lm_y_cathode_ceramic metal_interlayer_normal_lm metal_pore_normal_lm interlayer_pore_normal_lm"
  # group_variables = 'ux uy normal_l
  acceptable_iterations = 2
[]

[Outputs]
  exodus = true
  execute_on = 'INITIAL NONLINEAR TIMESTEP_END'
[]

[Debug]
  show_var_residual_norms = true
  # show_material_props = true
[]
