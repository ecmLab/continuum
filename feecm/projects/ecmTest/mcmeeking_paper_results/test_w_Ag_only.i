# Test to check results of mcmeeking paper
# Test overpotential on the surface etc.

current_density = 10e-3
k_anode = 1e-3
faraday = 96.4853329
temperature = 298
gas_constant = 8.31446
c_ref = 0
k_pore = 1e-3


c_init_Carbon = 2.5925e-2
c_max_Carbon = 3.05e-2
c_ref_Carbon = 2.5925e-2

c_init_Silver = 72.2e-3
c_max_Silver = 173.2e-3
c_ref_Silver = 72.2e-3

c_init_Li = 0.0
c_ref_Li = 0.0
c_max_Li = 1.0

pressure = 1e-4
end_time = 150
[Mesh]
  construct_node_list_from_side_list = true
  [./mesh]
    type = FileMeshGenerator
    file = 'data/test_w_Ag_interlayer_only1.e'
  [../]
  #-------------------------------------------------------------#
  # --- interface between interlayer and ceramic
  [./secondary_Ag_block]
    type = LowerDBlockFromSidesetGenerator
    input = mesh
    sidesets = 'blockSilver_bottom'
    new_block_name = 'Ag_s_block'
  [../]
  [./primary_Ag_block]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_Ag_block
    sidesets = 'blockCeramic_top'
    new_block_name = 'Ag_p_block'
  [../]
  # --- interface between interLayer and Ag
  [./interLayer_Ag_secondary]
    type = LowerDBlockFromSidesetGenerator
    input = primary_Ag_block
    sidesets = 'interLayer_bottom'
    new_block_name = 'interLayer_s_block'
  [../]
  [./interLayer_Ag_primary]
    type = LowerDBlockFromSidesetGenerator
    input = interLayer_Ag_secondary
    sidesets = 'blockSilver_top'
    new_block_name = 'interLayer_p_block'
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
    block = 'blockCeramic blockSilver interLayer'
    extra_vector_tags = 'ref'
  [../]
[]

[Functions]
  [./pressure]
    type = ParsedFunction
    # value = 'if (t<1, 1e-2*t, 1e-2)'
    # value = 'if(t<=10, 1e1*t, 1e2)'
    value = 'if (t < 100, ${pressure}, (t-100)*1e1)'
  [../]
  [./current_density]
    type = ParsedFunction
    # value = 'if (t > 10, ${current_density},1e-6)'
    value = 'if (t <=10, ${current_density}/10.0*t, ${current_density})'
  [../]
  [./reaction_rate_c]
    type = PiecewiseLinear
    data_file = 'data/carbon_black_equil_potential_mV.csv'
    format = columns
  [../]
[]



[Variables]
  [./ux]
    block = 'blockCeramic blockSilver interLayer'
  [../]
  [./uy]
    block = 'blockCeramic blockSilver interLayer'
  [../]
  [./V]
    block = 'blockCeramic blockSilver'
  [../]
  [./li_conc]
    block = 'blockSilver interLayer'
  [../]

  # Ion flux at ag ceramic boundary
  [./flux_Ag_ceramic]
    block = 'Ag_s_block'
  [../]
  [./normal_lm_x_Ag]
    block = 'Ag_s_block'
  [../]
  [./normal_lm_y_Ag]
    block = 'Ag_s_block'
  [../]
  # --- interlayer Flux
  [./flux_interlayer_Ag]
    block = 'interLayer_s_block'
  [../]
  [./normal_lm_x_interlayer]
    block = 'interLayer_s_block'
  [../]
  [./normal_lm_y_interLayer]
    block = 'interLayer_s_block'
  [../]

[]

[ICs]
  [./c_init_Ag]
    type = ConstantIC
    value = ${c_init_Silver}
    block = 'blockSilver'
    variable = li_conc
  [../]
  [./c_init_interlayer]
    type = ConstantIC
    value = ${c_init_Li}
    block = 'interLayer'
    variable = li_conc
  [../]
[]

[AuxVariables]
  [./chemical_potential]
    order = CONSTANT
    family = MONOMIAL
    block = 'blockSilver interLayer'
  [../]

  [./bndliflux]
    order = CONSTANT
    family = MONOMIAL
    block = 'blockCeramic blockSilver '
  [../]
  [./li_ion_flux_x]
    order = FIRST
    family = MONOMIAL
    block = 'blockCeramic blockSilver '
  [../]
  [./li_ion_flux_y]
    order = FIRST
    family = MONOMIAL
    block = 'blockCeramic blockSilver '
  [../]
  [./li_ion_flux_z]
    order = FIRST
    family = MONOMIAL
    block = 'blockCeramic blockSilver '
  [../]
  [./eq_pot]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./soc]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./li_metal_flux_x]
    order = CONSTANT
    family = MONOMIAL
    block = 'blockSilver interLayer'
  [../]
  [./li_metal_flux_y]
    order = CONSTANT
    family = MONOMIAL
    block = 'blockSilver interLayer'
  [../]
  [./li_metal_flux_z]
    order = CONSTANT
    family = MONOMIAL
    block = 'blockSilver interLayer'
  [../]
[]

[AuxKernels]
  [./li_metal_flux_x1]
    type = ADDiffusionFluxAux
    variable = li_metal_flux_x
    component = x
    diffusion_variable = li_conc
    block = 'blockSilver'
    diffusivity = diff_Ag
  [../]

  [./li_metal_flux_y1]
    type = ADDiffusionFluxAux
    variable = li_metal_flux_y
    component = y
    diffusion_variable = li_conc
    block = 'blockSilver'
    diffusivity = diff_Ag
  [../]
  [./li_metal_flux_z1]
    type = ADDiffusionFluxAux
    variable = li_metal_flux_z
    component = z
    diffusion_variable = li_conc
    block = 'blockSilver'
    diffusivity = diff_Ag
  [../]
  [./li_metal_flux_x2]
    type = AnisoTropicDiffusionFluxAux
    variable = li_metal_flux_x
    component = x
    diffusion_variable = li_conc
    block = 'interLayer'
    diffusivity = diffusivity
  [../]

  [./li_metal_flux_y2]
    type = AnisoTropicDiffusionFluxAux
    variable = li_metal_flux_y
    component = y
    diffusion_variable = li_conc
    block = 'interLayer'
    diffusivity = diffusivity
  [../]
  [./li_metal_flux_z2]
    type = AnisoTropicDiffusionFluxAux
    variable = li_metal_flux_z
    component = z
    diffusion_variable = li_conc
    block = 'interLayer'
    diffusivity = diffusivity
  [../]

  [./chemical_potential]
    type = ADMaterialRealAux
    variable = chemical_potential
    block = 'blockSilver interLayer'
    property = chemical_potential
  [../]
  [./bnd_li_flux]
    type = DiffusionFluxNormalToBoundaryAux
    variable = bndliflux
    boundary = 'blockCeramic_top blockSilver_bottom'
    diffusion_variable = V
    diffusivity = thermal_conductivity
  [../]
  [./li_ion_flux_x]
    type = ADDiffusionFluxAux
    variable = li_ion_flux_x
    component = x
    diffusion_variable = V
    diffusivity = thermal_conductivity
    block = 'blockCeramic blockSilver'
  [../]

  [./li_ion_flux_y]
    type = ADDiffusionFluxAux
    variable = li_ion_flux_y
    component = y
    diffusion_variable = V
    diffusivity = thermal_conductivity
    block = 'blockCeramic blockSilver '
  [../]
  [./li_ion_flux_z]
    type = ADDiffusionFluxAux
    variable = li_ion_flux_z
    component = z
    diffusion_variable = V
    diffusivity = thermal_conductivity
    block = 'blockCeramic blockSilver'
  [../]

  [./eq_pot]
    type = ADMaterialRealAux
    variable = eq_pot
    block = 'blockSilver'
    property = 'equilibrium_potential'
  [../]

  [./soc]
    type = ADMaterialRealAux
    variable = soc
    block = 'blockSilver'
    property = 'state_of_charge'
  [../]
[]

[Constraints]
  ### - Constraint between Ag Layer and ceramic
  ## -- All electrochemistry happens here
  [./AgLayer_constraint]
    type = GapDisplacementConductanceConstraint
    variable = flux_Ag_ceramic
    secondary_variable = V
    secondary_boundary =  'blockSilver_bottom'
    secondary_subdomain = 'Ag_s_block'
    primary_boundary = 'blockCeramic_top'
    primary_subdomain = 'Ag_p_block'
    concentration = 'li_conc'
    k = ${k_anode}
    use_displaced_mesh = true
    displacements = 'ux uy'
    compute_lm_residuals = true
    include_equilibrium_potential = false
    include_gap = false
    include_concentration = true
    faraday = ${faraday}
    temperature = ${temperature}
    extra_vector_tags = 'ref'
    surfaceType = SECONDARY
  [../]
  [./conc_Aglayer]
    type = ScaledBCConstraint
    variable = flux_Ag_ceramic
    secondary_variable = li_conc
    primary_variable = V
    primary = false
    secondary_boundary =  'blockSilver_bottom'
    secondary_subdomain = 'Ag_s_block'
    primary_boundary = 'blockCeramic_top'
    primary_subdomain = 'Ag_p_block'
    scale = 1.036426e-2
    use_displaced_mesh = true
    extra_vector_tags = 'ref'
  [../]
  [./Ag_mech_x]
    type = EqualValueConstraint
    variable = normal_lm_x_Ag
    primary_variable = ux
    secondary_variable = ux
    secondary_boundary =  'blockSilver_bottom'
    secondary_subdomain = 'Ag_s_block'
    primary_boundary = 'blockCeramic_top'
    primary_subdomain = 'Ag_p_block'
    use_displaced_mesh = true
  [../]
  [./Ag_mech_y]
    type = EqualValueConstraint
    variable = normal_lm_y_Ag
    primary_variable = uy
    secondary_variable = uy
    secondary_boundary =  'blockSilver_bottom'
    secondary_subdomain = 'Ag_s_block'
    primary_boundary = 'blockCeramic_top'
    primary_subdomain = 'Ag_p_block'
    use_displaced_mesh = true
  [../]
  # --- -Constraints between interLayer and Ag layer
  [./interlayer_Ag_mech_x]
    type = EqualValueConstraint
    variable = normal_lm_x_interlayer
    primary_variable = ux
    secondary_variable = ux
    secondary_boundary =  'interLayer_bottom'
    secondary_subdomain = 'interLayer_s_block'
    primary_boundary = 'blockSilver_top'
    primary_subdomain = 'interLayer_p_block'
    use_displaced_mesh = true
  [../]
  [./interlayer_Ag_mech_y]
    type = EqualValueConstraint
    variable = normal_lm_y_interLayer
    primary_variable = uy
    secondary_variable = uy
    secondary_boundary =  'interLayer_bottom'
    secondary_subdomain = 'interLayer_s_block'
    primary_boundary = 'blockSilver_top'
    primary_subdomain = 'interLayer_p_block'
    use_displaced_mesh = true
  [../]
  # [./equilibriate_interlayer]
  #   type = EqualGradientConstraint
  #   variable = flux_interlayer_Ag
  #   secondary_boundary =  'interLayer_bottom'
  #   secondary_subdomain = 'interLayer_s_block'
  #   primary_boundary = 'blockSilver_top'
  #   primary_subdomain = 'interLayer_p_block'
  #   primary_variable = 'li_conc'
  #   secondary_variable = 'li_conc'
  #   component = 1
  #   extra_vector_tags = 'ref'
  # [../]
  # [./equilibriate_interlayer]
  #   type = GapEquilibriateConstraint
  #   variable = flux_interlayer_Ag
  #   primary_mat_prop = chemical_potential
  #   secondary_mat_prop = chemical_potential
  #   primary_variable = li_conc
  #   secondary_variable = li_conc
  #   k = -1e-6
  #   use_displaced_mesh = true
  #   include_gap = false
  #   temperature = ${temperature}
  #   faraday = ${faraday}
  #   R = ${gas_constant}
  #   secondary_boundary =  'interLayer_bottom'
  #   secondary_subdomain = 'interLayer_s_block'
  #   primary_boundary = 'blockSilver_top'
  #   primary_subdomain = 'interLayer_p_block'
  #   extra_vector_tags = 'ref'
  #   prefactor = RT
  #   # one_sided = SECONDARY->PRIMARY
  # [../]
  [./equilibriate_interlayer]
    type = EqualNormalFluxConstraint
    secondary_boundary =  'interLayer_bottom'
    secondary_subdomain = 'interLayer_s_block'
    primary_boundary = 'blockSilver_top'
    primary_subdomain = 'interLayer_p_block'
    extra_vector_tags = 'ref'
    primary_variable = 'li_conc'
    secondary_variable = 'li_conc'
    primary_mat_prop = diff_Ag
    secondary_mat_prop = diffusivity
    primary_tensor = false
    secondary_tensor = true
    variable = flux_interlayer_Ag
  [../]
[]

[Kernels]
  [./V]
    type = ADHeatConduction
    variable = V
    use_displaced_mesh = false
    block = 'blockCeramic blockSilver'
  [../]
  [./li_metal_interlayer]
    type = ADChemoMechanoAnsioDiffusion
    variable = li_conc
    diffusivity = diffusivity
    use_displaced_mesh = false
    block = 'interLayer'
  [../]
  [./li_metal_Ag]
    type = ADChemoMechanoDiffusion
    variable = li_conc
    diffusivity = diff_Ag
    use_displaced_mesh = false
    block = 'blockSilver'
  [../]
  [./li_metal_dt]
    type = ADTimeDerivative
    variable = li_conc
    use_displaced_mesh = false
    block = 'blockSilver interLayer'
  [../]
[]

[Materials]
  [./thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 100
    block = 'blockSilver interLayer'
  [../]
  [./thermal_conductivity2]
    type = ADHeatConductionMaterial
    thermal_conductivity = 0.1
    block = 'blockCeramic'
  [../]
  [./diff_in_interlayer]
    type = ADDiffusionAlongPrincipalDirectionsMaterial
    diffusivity_vector = '5e6 0 0'
    block = 'interLayer'
  [../]
  [./diff_in_Ag]
    type = ADGenericConstantMaterial
    prop_names = 'diff_Ag'
    prop_values = '1000'
    block = 'blockSilver'
  [../]
  [./equilibrium_potential_Li]
    type = ADComputeEquilibriumPotential
    R = ${gas_constant}
    faraday = ${faraday}
    temperature = ${temperature}
    cref = ${c_max_Li}
    concentration = li_conc
    include_conc = false
    include_reaction_rate = true
    reaction_rate = 0.0
    # reaction_rate_function = reaction_rate
    include_mechanical_effects = false
    exclude_elastic_contribution = true
    block = 'interLayer'
  [../]
  [./equilibrium_potential_silver]
    type = ADComputeEquilibriumPotential
    R = ${gas_constant}
    faraday = ${faraday}
    temperature = ${temperature}
    cref = ${c_max_Silver}
    concentration = li_conc
    include_conc = true
    include_reaction_rate = true
    reaction_rate = 0
    # reaction_rate_function = reaction_rate_c
    include_mechanical_effects = false
    exclude_elastic_contribution = true
    block = 'blockSilver'
    potential = BOWER
  [../]
  # ----- Silver Mechanical Properties
  [./elasticity_tensor_Ag]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 7.8e3
    poissons_ratio = 0.3
    block = 'blockSilver'
  [../]
  [./elasticity_tensor_Llzo]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 170e3
    poissons_ratio = 0.25
    block = 'blockCeramic'
  [../]
  [./stress_llzo]
    type = ADComputeFiniteStrainElasticStress
    block = 'blockCeramic'
  [../]
  [./stress_Li_Ag]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas_Ag'
    # perform_finite_strain_rotations = true
    block = 'blockSilver'
  [../]

  [./plas_Ag]
    type = ADIsoTropicHyperViscoSwellingCreep
    rate_form = false
    absolute_tolerance = 1e-8
    relative_tolerance = 1e-6
    isotropic_swelling = true
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
    alpha = '0.33333 0.33333 0.33333'
    concentration = li_conc
    cref = ${c_ref_Silver}
    block = 'blockSilver'
    intBnd = 'blockCeramic_top'
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    enable = true
  [../]
  # Interlayer Properties
  [./elasticity_tensor_Li]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 7.8e3
    poissons_ratio = 0.3
    block = 'interLayer'
  [../]

  [./stress_Li]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas_Li'
    # perform_finite_strain_rotations = true
    block = 'interLayer'
  [../]
  [./plas_Li]
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
    concentration = li_conc
    cref = ${c_ref_Li}
    block = 'interLayer'
    intBnd = 'blockCeramic_top'
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    enable = true
  [../]
[]

[BCs]
  [./current]
    type = ADFunctionNeumannBC
    boundary = 'blockCeramic_bottom'
    # value = ${current_density}
    function = current_density
    variable = V
    extra_vector_tags = 'ref'
  [../]

  [./OV]
    type = ADDirichletBC
    boundary = 'blockSilver_top'
    value = 0.0
    variable = V
  [../]

  [./ux]
    type = ADDirichletBC
    preset = true
    boundary = 'blockCeramic_left blockCeramic_right blockSilver_left blockSilver_right interLayer_left interLayer_right'
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
    boundary = 'interLayer_top'
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
  dt = 200
  # num_steps = 20
  l_max_its = 50
  nl_max_its = 25
  nl_abs_tol = 1e-11
  nl_rel_tol = 1e-3
  dtmax = 200
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 15
    # timestep_limiting_postprocessor = matl_ts_min
  [../]
  # num_steps = 1
  end_time = ${end_time}
  scaling_group_variables = 'V flux_Ag_ceramic'
  resid_vs_jac_scaling_param = 0.5
[]
[Problem]
  type = ReferenceResidualProblem
  # solution_variables = 'ux uy normal_lm thermal_lm V thermal_lm2 li_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'V flux_Ag_ceramic li_conc flux_interlayer_Ag; ux uy normal_lm_x_Ag normal_lm_y_Ag normal_lm_x_interlayer normal_lm_y_interlayer'
  # group_variables = 'V flux_interlayer_Ag flux_pore_Ag li_conc; ux uy normal_lm_x_Ag normal_lm_y_Ag normal_lm_x_interlayer normal_lm_y_interlayer'
  acceptable_iterations = 2
  # coord_type = RZ
[]
[Outputs]
  # csv = true
  [./csv]
    type = CSV
    file_base = csv/test_pore
  [../]
  [./out]
    type = Exodus
    file_base = rst/test_pore
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
[]

[Postprocessors]
  [./chemical_pot_interlayer]
    type = SideAverageValue
    boundary = 'interLayer_bottom'
    variable = chemical_potential
  [../]
  [./chemical_pot_Ag]
    type = SideAverageValue
    boundary = 'blockSilver_top'
    variable = chemical_potential
  [../]

  [./bot_conc_Li]
    type = SideAverageValue
    boundary = 'interLayer_bottom'
    variable = li_conc
  [../]
  [./top_conc_silver]
    type = SideAverageValue
    boundary = 'blockSilver_top'
    variable = li_conc
  [../]
  [./flux_interlayer_Ag]
    type = AverageNodalVariableValue
    variable = flux_interlayer_Ag
    block = 'interLayer_s_block'
  [../]

[]

[VectorPostprocessors]
  [./interLayer_Ag_flux]
    type = NodalValueSampler
    sort_by = x
    block = 'interLayer_s_block'
    variable = flux_interlayer_Ag
  [../]
[]
