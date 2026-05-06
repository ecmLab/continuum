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
pressure = 1.0
end_time = 100.0
[Mesh]
  construct_node_list_from_side_list = true
  patch_update_strategy = iteration
  patch_size = 100
  # file = test_pore_equil_cp/0025_mesh.cpr
  [./mesh]
    type = FileMeshGenerator
    # file = 'data/interlayer_ceramic_metal_pore2.e'
    file = 'data/test_w_pore_interlayer_w_block_metal_h0.5_shallow.e'
    # file = 'data/interlayer_w_top.e'
    # file = 'data/test_pore5.e'
  [../]
  # --- interface between pore and ceramic
  [./pore_ceramic_secondary]
    type = LowerDBlockFromSidesetGenerator
    input = mesh
    sidesets = 'blockPore_bottom'
    new_block_name = 'pore_s_block'
  [../]
  [./pore_ceramic_primary]
    type = LowerDBlockFromSidesetGenerator
    input = pore_ceramic_secondary
    sidesets = 'blockCeramic_top'
    new_block_name = 'pore_p_block'
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
    block = 'blockCeramic interLayer blockMetal blockPore'
    extra_vector_tags = 'ref'
  [../]
[]

[Contact]
  [./interlayer_ceramic]
    formulation = mortar
    model = frictionless
    primary =  'blockCeramic_top'
    secondary =  'interLayer_bottom'
    normal_smoothing_distance = 0.5
    tangential_tolerance = 0.25
    mesh = pore_ceramic_primary
  [../]

  [./metal_interlayer]
    formulation = mortar
    model = frictionless
    primary =  'blockMetal_bottom'
    secondary =  'interLayer_top'
    normal_smoothing_distance = 0.5
    tangential_tolerance = 0.25
    mesh = interlayer_ceramic_secondary_subdomain_generator
  [../]
  [./metal_pore]
    formulation = mortar
    model = frictionless
    primary =  'blockPore_left'
    secondary =  'blockMetal_right'
    normal_smoothing_distance = 0.5
    tangential_tolerance = 0.25
    mesh = metal_interlayer_secondary_subdomain_generator
  [../]
  [./interlayer_pore]
    formulation = mortar
    model = frictionless
    primary =  'blockPore_left'
    secondary =  'interLayer_right'
    normal_smoothing_distance = 0.5
    tangential_tolerance = 0.25
    mesh = metal_pore_secondary_subdomain_generator
  [../]
  # [./pore_ceramic]
  #   formulation = mortar
  #   model = frictionless
  #   primary =  'blockCeramic_top'
  #   secondary =  'blockPore_bottom'
  #   normal_smoothing_distance = 0.5
  #   tangential_tolerance = 0.25
  #   mesh = interlayer_pore_secondary_subdomain_generator
  # [../]

[]

[Functions]
  [./pressure]
    type = ParsedFunction
    # value = 'if (t<1, 1e-2*t, 1e-2)'
    # value = 'if(t<=10, 1e1*t, 1e2)'
    value = 'if (t <= 100, ${pressure}*t, if (t<=110,-9.95*(t-100)+100, 0.5))'
    # value = '${pressure}'
  [../]
  [./current_density]
    type = ParsedFunction
    # value = 'if (t > 10, ${current_density},1e-6)'
    value = 'if (t <=10, if (t < 1, 0, ${current_density}/10.0*t), ${current_density})'
  [../]
  [./reaction_rate_c]
    type = PiecewiseLinear
    data_file = 'data/carbon_black_equil_potential_mV.csv'
    format = columns
  [../]
[]

[Variables]
  [./ux]
    block = 'blockCeramic interLayer blockMetal blockPore'
  [../]
  [./uy]
    block = 'blockCeramic interLayer blockMetal blockPore'
  [../]
  [./normal_lm_x_pore]
    block = 'pore_s_block'
  [../]
  [./normal_lm_y_pore]
    block = 'pore_s_block'
  [../]
  # [./normal_lm_x_anode]
  #     block = 'anode_s_block'
  # [../]
  # [./normal_lm_y_anode]
  #   block = 'anode_s_block'
  # [../]
[]


[AuxVariables]
  [./plastic_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./li_conc2]
    block = 'interLayer'
  [../]
[]

[AuxKernels]
  [./li_conc2]
    type = ConstantAux
    value = 0.0
    variable = li_conc2
    block = 'interLayer'
  [../]

  [./peeq]
    type = ADMaterialRealAux
    variable = plastic_strain
    property = effective_plastic_strain
    execute_on = timestep_end
    block = 'interLayer'
  [../]
[]

[Constraints]
  [./pore_mech_x]
    type = EqualValueConstraint
    variable = normal_lm_x_pore
    primary_variable = ux
    secondary_variable = ux
    secondary_boundary =  'blockPore_bottom'
    secondary_subdomain = 'pore_s_block'
    primary_boundary = 'blockCeramic_top'
    primary_subdomain = 'pore_p_block'
    use_displaced_mesh = true
  [../]
  [./pore_mech_y]
    type = EqualValueConstraint
    variable = normal_lm_y_pore
    primary_variable = uy
    secondary_variable = uy
    secondary_boundary =  'blockPore_bottom'
    secondary_subdomain = 'pore_s_block'
    primary_boundary = 'blockCeramic_top'
    primary_subdomain = 'pore_p_block'
    use_displaced_mesh = true
  [../]
[]


[Materials]
  [./equilibrium_potential_Li]
    type = ADComputeEquilibriumPotential
    R = ${gas_constant}
    faraday = ${faraday}
    temperature = ${temperature}
    cref = 1.0
    concentration = li_conc2
    include_conc = false
    include_reaction_rate = true
    reaction_rate = 0.0
    # reaction_rate_function = reaction_rate
    include_mechanical_effects = false
    exclude_elastic_contribution = true
    block = 'interLayer'
  [../]

  [./elasticity_tensor_Li]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 7.8e3
    poissons_ratio = 0.3
    block = 'interLayer blockMetal'
  [../]
  [./elasticity_tensor_C]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 100e3
    poissons_ratio = 0.3
    block = 'blockPore'
  [../]
  [./elasticity_tensor_Llzo]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 170e3
    poissons_ratio = 0.25
    block = 'blockCeramic'
  [../]
  [./stress_llzo]
    type = ADComputeFiniteStrainElasticStress
    block = 'blockCeramic blockPore'
  [../]
  [./stress_Li]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    # perform_finite_strain_rotations = true
    block = 'interLayer'
  [../]
  [./plas]
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
    omega = 13
    alpha = '1 0 0'
    concentration = li_conc2
    cref = ${c_ref}
    block = 'interLayer'
    intBnd = 'blockCeramic_top'
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    enable = true
  [../]
  [./stress_Li2]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas2'
    # perform_finite_strain_rotations = true
    block = 'blockMetal'
  [../]
  [./plas2]
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
  [../]

[]

[BCs]
  [./ux]
    type = ADDirichletBC
    preset = true
    boundary = 'blockCeramic_left blockCeramic_right interLayer_left blockMetal_left blockPore_right'
    variable = ux
    value = 0
  [../]

  [./uy]
    type = ADDirichletBC
    preset = true
    boundary = 'blockCeramic_bottom'# blockMetal_top'
    variable = uy
    value = 0
  [../]
  [./pressure]
    type = ADPressure
    variable = uy
    boundary = 'blockMetal_top'
    component = 1
    function = pressure
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
  nl_max_its = 25
  nl_abs_tol = 1e-11
  nl_rel_tol = 1e-3
  dtmax = 5
  dtmin = 1e-6
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-2
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 15
    timestep_limiting_postprocessor = matl_ts_min
  [../]
  # num_steps = 1
  end_time = ${end_time}
  # scaling_group_variables = 'V flux_interlayer'# flux_pore'
  resid_vs_jac_scaling_param = 0.5
[]
[Problem]
  type = ReferenceResidualProblem
  # solution_variables = 'ux uy normal_lm thermal_lm V thermal_lm2 li_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'ux uy  metal_interlayer_normal_lm metal_pore_normal_lm interlayer_pore_normal_lm normal_lm_x_pore normal_lm_y_pore interlayer_ceramic_normal_lm'
  acceptable_iterations = 2
  # restart_file_base = test_pore_equil_cp/0025
  # coord_type = RZ
[]
[Outputs]
  # csv = true
  # exodus = true
  # execute_on = 'NONLINEAR'
  [./csv]
    type = CSV
    file_base = csv/pore
  [../]
  [./out]
    type = Exodus
    file_base = rst/pore
    execute_on = 'INITIAL TIMESTEP_END'
    # interval = 1
  [../]
  [./check]
    type = Checkpoint
    file_base = check/pore
    start_time = 1.0
    interval = 5
    num_files = 4
  [../]
[]

[Postprocessors]
  # [./eq_pot]
  #   type = SideAverageValue
  #   boundary = 'interLayer_bottom'
  #   variable = eq_pot
  # [../]
  # [./eq_pot_pore]
  #   type = SideAverageValue
  #   boundary = 'blockPore_bottom'
  #   variable = eq_pot
  # [../]
  [./bot_conc_Li]
    type = SideAverageValue
    boundary = 'interLayer_bottom'
    variable = li_conc2
  [../]
  # [./bot_conc_pore]
  #   type = SideAverageValue
  #   boundary = 'blockPore_bottom'
  #   variable = li_conc
  # [../]

  [./matl_ts_min]
    type = MaterialTimeStepPostprocessor
    block = 'interLayer blockMetal'
  [../]
  # [./pore_bot_soc]
  #   type = SideAverageValue
  #   variable = soc
  #   boundary = 'blockPore_bottom'
  # [../]
  # [./interlayer_pore_flux]
  #   type = ElementAverageValue
  #   variable = equil
  #   block = 'interlayer_pore_secondary_subdomain'
  # [../]
  # [./contact_size]
  #   type = ContactDOFSetSize
  #   variable = interlayer_ceramic_normal_lm
  #   subdomain = interlayer_ceramic_secondary_subdomain
  #   tolerance = 1e-10
  # [../]
  # [./contact_size_pore]
  #   type = ContactDOFSetSize
  #   variable = interlayer_pore_normal_lm
  #   subdomain = interlayer_pore_secondary_subdomain
  #   tolerance = 1e-10
  # [../]
[]
