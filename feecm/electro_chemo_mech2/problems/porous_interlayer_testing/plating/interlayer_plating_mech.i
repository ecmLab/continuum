current_density = 10e-3
k_anode = 1e-3
faraday = 96.4853329
temperature = 298
gas_constant = 8.31446
c_ref = 0
# c_ref_interLayer = 2.592e-2
c_ref_interLayer = 0.0
k_interLayer = 1e-3
c_init_interLayer = 0.0
# c_init_C = 2.592e-2
c_init_blockPore = 0.0
c_max_C = 3.05e-2
c_ref_blockPore = 0.0
pressure = 0.1
end_time = 1
metal_flux = 2.0534e-3
base_name = 'interlayer_mech'
[Mesh]
  construct_node_list_from_side_list = true
  patch_update_strategy = iteration
  patch_size = 100
  # file = check/interLayer_cp/0002_mesh.cpr
  [./mesh]
    type = FileMeshGenerator
    # file = 'data/interlayer_left_pore_right.e'
    file = 'data/test_w_pore_interlayer_w_block_metal_h0.5_shallow.e'
  [../]
  [./interLayer_ceramic_secondary]
    type = LowerDBlockFromSidesetGenerator
    input = mesh
    sidesets = 'interLayer_bottom'
    new_block_name = 'interLayer_s_block'
  [../]
  [./interLayer_ceramic_primary]
    type = LowerDBlockFromSidesetGenerator
    input = interLayer_ceramic_secondary
    sidesets = 'blockCeramic_top'
    new_block_name = 'interLayer_p_block'
  [../]

  [./blockPore_ceramic_secondary]
    type = LowerDBlockFromSidesetGenerator
    input = interLayer_ceramic_primary
    sidesets = 'blockPore_bottom'
    new_block_name = 'blockPore_s_block'
  [../]
  [./blockPore_ceramic_primary]
    type = LowerDBlockFromSidesetGenerator
    input = blockPore_ceramic_secondary
    sidesets = 'blockCeramic_top'
    new_block_name = 'blockPore_p_block'
  [../]

  [./blockMetal_interLayer_secondary]
    type = LowerDBlockFromSidesetGenerator
    input = blockPore_ceramic_primary
    sidesets = 'interLayer_top'
    new_block_name = 'blockMetal_interLayer_s_block'
  [../]
  [./blockMetal_interLayer_primary]
    type = LowerDBlockFromSidesetGenerator
    input = blockMetal_interLayer_secondary
    sidesets = 'blockMetal_bottom'
    new_block_name = 'blockMetal_interLayer_p_block'
  [../]
[]
[GlobalParams]
  displacements = 'ux uy'
[]

[Contact]
  [./interLayer_blockPore]
    formulation = mortar
    model = frictionless
    primary =  'blockPore_left'
    secondary =  'interLayer_right'
    normal_smoothing_distance = 0.5
    tangential_tolerance = 0.25
    mesh = blockMetal_interLayer_primary
    # mesh = mesh
  [../]
  [./blockMetal_blockPore]
    formulation = mortar
    model = frictionless
    primary =  'blockPore_left'
    secondary =  'blockMetal_right'
    normal_smoothing_distance = 0.5
    tangential_tolerance = 0.25
    mesh = interLayer_blockPore_secondary_subdomain_generator
    # mesh = mesh
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress hydrostatic_stress'
    use_automatic_differentiation = true
    block = 'blockCeramic interLayer blockPore blockMetal'
    extra_vector_tags = 'ref'
  [../]
[]

[Variables]
  [./ux]
    block = 'blockCeramic interLayer  blockMetal blockPore'
  [../]
  [./uy]
    block = 'blockCeramic  interLayer  blockMetal blockPore'
  [../]

  [./normal_lm_x_interLayer]
    block = 'interLayer_s_block'
  [../]
  [./normal_lm_y_interLayer]
    block = 'interLayer_s_block'
  [../]
  # ------------- Pore variables
  [./normal_lm_x_blockPore]
    block = 'blockPore_s_block'
  [../]
  [./normal_lm_y_blockPore]
    block = 'blockPore_s_block'
  [../]
  # ------------- interLayer metal variables
  [./normal_lm_x_blockMetal_interLayer]
    block = 'blockMetal_interLayer_s_block'
  [../]
  [./normal_lm_y_blockMetal_interLayer]
    block = 'blockMetal_interLayer_s_block'
  [../]
[]

[Constraints]
  [./interLayer_mech_x]
    type = EqualValueConstraint
    variable = normal_lm_x_interLayer
    primary_variable = ux
    secondary_variable = ux
    secondary_boundary =  'interLayer_bottom'
    secondary_subdomain = 'interLayer_s_block'
    primary_boundary = 'blockCeramic_top'
    primary_subdomain = 'interLayer_p_block'
    use_displaced_mesh = true
  [../]
  [./interLayer_mech_y]
    type = EqualValueConstraint
    variable = normal_lm_y_interLayer
    primary_variable = uy
    secondary_variable = uy
    secondary_boundary =  'interLayer_bottom'
    secondary_subdomain = 'interLayer_s_block'
    primary_boundary = 'blockCeramic_top'
    primary_subdomain = 'interLayer_p_block'
    use_displaced_mesh = true
  [../]
  # ---- Pore Ceramic glued constraint
  [./blockPore_mech_x]
    type = EqualValueConstraint
    variable = normal_lm_x_blockPore
    primary_variable = ux
    secondary_variable = ux
    secondary_boundary =  'blockPore_bottom'
    secondary_subdomain = 'blockPore_s_block'
    primary_boundary = 'blockCeramic_top'
    primary_subdomain = 'blockPore_p_block'
    use_displaced_mesh = true
  [../]
  [./blockPore_mech_y]
    type = EqualValueConstraint
    variable = normal_lm_y_blockPore
    primary_variable = uy
    secondary_variable = uy
    secondary_boundary =  'blockPore_bottom'
    secondary_subdomain = 'blockPore_s_block'
    primary_boundary = 'blockCeramic_top'
    primary_subdomain = 'blockPore_p_block'
    use_displaced_mesh = true
  [../]
  # ---- Metal interlayer glued constraint
  [./blockMetal_interLayer_mech_x]
    type = EqualValueConstraint
    variable = normal_lm_x_blockMetal_interLayer
    primary_variable = ux
    secondary_variable = ux
    secondary_boundary =  'interLayer_top'
    secondary_subdomain = 'blockMetal_interLayer_s_block'
    primary_boundary = 'blockMetal_bottom'
    primary_subdomain = 'blockMetal_interLayer_p_block'
    use_displaced_mesh = true
  [../]
  [./blockMetal_interLayer_mech_y]
    type = EqualValueConstraint
    variable = normal_lm_y_blockMetal_interLayer
    primary_variable = uy
    secondary_variable = uy
    secondary_boundary =  'interLayer_top'
    secondary_subdomain = 'blockMetal_interLayer_s_block'
    primary_boundary = 'blockMetal_bottom'
    primary_subdomain = 'blockMetal_interLayer_p_block'
    use_displaced_mesh = true
  [../]


[]

[Materials]
  [./thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 100
    block = 'interLayer blockPore'
  [../]
  [./thermal_conductivity2]
    type = ADHeatConductionMaterial
    thermal_conductivity = 0.1
    block = 'blockCeramic'
  [../]

  [./diff_in_interLayer]
    type = ADDiffusionAlongPrincipalDirectionsMaterial
    diffusivity_vector = '5e6 0.5 0.5'
    block = 'interLayer'
  [../]
  [./diff_in_blockPore]
    type = ADGenericConstantMaterial
    prop_names = 'diff_blockPore'
    prop_values = '0.15'
    block = 'blockPore'
  [../]
  [./equilibrium_potential_interlayer]
    type = ADComputeEquilibriumPotential
    R = ${gas_constant}
    faraday = ${faraday}
    temperature = ${temperature}
    cref = 0.1
    concentration = li_conc2
    include_conc = false
    include_reaction_rate = true
    # reaction_rate = 0.0
    is_ocv = true
    reaction_rate_function = reaction_rate_ag
    include_mechanical_effects = true
    exclude_elastic_contribution = true
    block = 'interLayer'
  [../]


  [./equilibrium_potential_blockPore]
    type = ADComputeEquilibriumPotential
    R = ${gas_constant}
    faraday = ${faraday}
    temperature = ${temperature}
    cref = ${c_max_C}
    concentration = li_conc2
    include_conc = false
    include_reaction_rate = true
    # reaction_rate = 0.0
    is_ocv = true
    reaction_rate_function = reaction_rate_c
    include_mechanical_effects = true
    exclude_elastic_contribution = true
    block = 'blockPore'
    # potential = DEHOFF
  [../]

  [./elasticity_tensor_Ceramic]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 170e3
    poissons_ratio = 0.25
    block = 'blockCeramic'
  [../]
  [./stress_Ceramic]
    type = ADComputeFiniteStrainElasticStress
    block = 'blockCeramic'
  [../]

  [./elasticity_tensor_interLayer]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 7.8e3
    poissons_ratio = 0.3
    block = 'interLayer blockMetal'
  [../]
  [./stress_interLayer]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plasticity_interLayer'
    # perform_finite_strain_rotations = true
    block = 'interLayer'
  [../]
  [./plasticity_interLayer]
    type = ADIsoTropicHyperViscoSwellingCreep
    rate_form = false
    absolute_tolerance = 1e-8
    relative_tolerance = 1e-6
    hardening_exponent = 2.0
    saturation_resistance = 2.0
    initial_resistance = 0.95
    hardening_modulus = 10.0
    rate_exponent = 0.15
    # reference_strain_rate = 0.05
    activation_energy = 37e3
    gas_constant = ${gas_constant} # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = ${temperature}

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 2.0
    omega = 1e-8
    alpha = '1 0 0'
    concentration = li_conc2
    cref = 0.0
    block = 'interLayer'
    intBnd = 'blockCeramic_top'
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    enable = true
  [../]
  [./elasticity_tensor_blockPore]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 20e3
    poissons_ratio = 0.3
    block = 'blockPore'
  [../]

  [./stress_blockPore]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'elastic_blockPore'
    # perform_finite_strain_rotations = true
    block = 'blockPore'
  [../]
  [./elastic_blockPore]
    type = ADIsotropicElasticSwelling
    # rate_form = false
    isotropic_swelling = true
    absolute_tolerance = 1e-8
    relative_tolerance = 1e-6
    alpha = '0.333333 0.333333 0.333333'
    omega = 1e-8
    block = 'blockPore'
    enable = true
    cref = ${c_ref_blockPore}
    concentration = li_conc2
  [../]
  [./stress_blockMetal]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plasticity_blockMetal'
    # perform_finite_strain_rotations = true
    block = 'blockMetal'
  [../]
  [./plasticity_blockMetal]
    type = ADIsoTropicHyperViscoCreep
    rate_form = false
    absolute_tolerance = 1e-8
    relative_tolerance = 1e-6
    hardening_exponent = 2.0
    saturation_resistance = 2.0
    initial_resistance = 0.95
    hardening_modulus = 10.0
    rate_exponent = 0.15
    # reference_strain_rate = 0.05
    activation_energy = 37e3
    gas_constant = ${gas_constant} # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = ${temperature}

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 2.0
    block = 'blockMetal'
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    enable = true
  [../]
[]

[BCs]

  [./ux]
    type = ADDirichletBC
    preset = true
    boundary = 'blockCeramic_left blockCeramic_right interLayer_left blockPore_right blockMetal_left'
    variable = ux
    value = 0
  [../]

  [./uy]
    type = ADDirichletBC
    preset = true
    boundary = 'blockCeramic_bottom'
    variable = uy
    value = 0
  [../]
  [./pressure]
    type = ADPressure
    variable = uy
    boundary = 'blockMetal_top blockPore_top'
    component = 1
    function = pressure
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
  solve_type = 'NEWTON'
  automatic_scaling = true
  # compute_scaling_once = false
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1     NONZERO               1e-20               '
  dt = 50
  # num_steps = 20
  l_max_its = 50
  nl_max_its = 15
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-5
  dtmax = 5
  dtmin = 1e-6
  nl_forced_its = 2
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.2
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 15
    timestep_limiting_postprocessor = matl_ts_min
  [../]
  # num_steps = 1
  end_time = ${end_time}
  # scaling_group_variables = 'V flux_interlayer flux_blockPore'
  resid_vs_jac_scaling_param = 0.5
[]
[Problem]
  type = ReferenceResidualProblem
  # solution_variables = 'ux uy normal_lm thermal_lm V thermal_lm2 li_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'ux uy normal_lm_x_interLayer normal_lm_y_interLayer normal_lm_x_blockPore normal_lm_y_blockPore interLayer_blockPore_normal_lm blockMetal_interLayer_normal_lm blockMetal_blockPore_normal_lm'#; V flux_interlayer flux_blockPore li_conc'
  # group_variables = 'ux uy pore_ceramic_normal_lm interLayer_ceramic_normal_lm interLayer_blockPore_normal_lm'#; V flux_interlayer flux_blockPore li_conc'
  acceptable_iterations = 2
  # restart_file_base = check/interLayer_cp/0002
  # coord_type = RZ
[]


[Outputs]
  # csv = true
  exodus = true
  execute_on = 'NONLINEAR'
  [./csv]
    type = CSV
    file_base = csv/${base_name}
  [../]
  [./out]
    type = Exodus
    file_base = rst/${base_name}
    execute_on = 'INITIAL TIMESTEP_END'
    sync_times = '1.0 10.0'
    sync_only = false
    # interval = 1
  [../]
  [./check]
    type = Checkpoint
    file_base = check/${base_name}
    start_time = 0.0
    interval = 5
    num_files = 4
    sync_times = '1.0 10.0'
    sync_only = false
  [../]
[]

[Functions]
  [./pressure]
    type = ParsedFunction
    # value = 'if (t<1, 1e-2*t, 1e-2)'
    # value = 'if(t<=10, 1e1*t, 1e2)'
    value = 'if (t <= 1, ${pressure}*t, ${pressure})'
    # value = '${pressure}'
  [../]
  [./current_density]
    type = ParsedFunction
    # value = 'if (t > 10, ${current_density},1e-6)'
    value = 'if (t <=10, if (t <= 1, 0.0, ${current_density}/9.0*(t-1)), ${current_density})'
  [../]
  [./reaction_rate_c]
    type = PiecewiseLinear
    data_file = 'data/carbon_black_equil_potential_mV.csv'
    format = columns
  [../]
  [./reaction_rate_ag]
    type = ParsedFunction
    value = 'if (t > 1.0, 0, 350.0 -350.0*t)'
  [../]
  [./interlayer_conc]
    type = ParsedFunction
    # value = 'if (t > 10, 1e-3*(t-10), 0)'
    value = 'if (t <=1, 0, if (t > 10, 5e-2, 5e-2/9 * (t-1) ))'
  [../]
  [./pore_conc]
    type = ParsedFunction
    # value = 'if (t > 1, 2.5925e-2/9 * (t-1), if (t > 10, 2.5925e-2, 0))'
    value = 'if (t <=1, 0, if (t > 10, 2.5925e-2, 2.5925e-2/9 * (t-1) ))'
  [../]
[]
[AuxVariables]
  [./li_conc2]
    block = 'interLayer blockPore'
  [../]
[]

[AuxKernels]
  [./li_conc2]
    type = FunctionAux
    function = interlayer_conc
    # value = 0
    variable = li_conc2
    block = 'interLayer'
  [../]

  [./li_conc]
    type = FunctionAux
    function = pore_conc
    variable = li_conc2
    block = 'blockPore'
  [../]

[]

[Postprocessors]
  [./matl_ts_min]
    type = MaterialTimeStepPostprocessor
    block = 'interLayer'
  [../]
[]
