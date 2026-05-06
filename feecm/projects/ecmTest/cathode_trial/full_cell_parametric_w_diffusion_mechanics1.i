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

current_density = 20e-3
k_anode = 1.0e-04
faraday = 96.4853329
temperature = 333
gas_constant = 8.31446
c_ref_cat = 4.9e-2
c_ref_anode = 0.1
c_init_cat = 2.45e-3
c_init_anode = 0.1
v_cutoff = 2750
h = 26

[Mesh]
  [./mesh]
    type = FileMeshGenerator
    file = 'test_full4.msh'
  [../]
  [./primary_boundary_cathode]
    type = LowerDBlockFromSidesetGenerator
    input = mesh
    sidesets = 'Electrolyte_bottom'
    new_block_name = 'cathode_p_block'
  [../]
  [./secondary_boundary_cathode]
    type = LowerDBlockFromSidesetGenerator
    input = primary_boundary_cathode
    sidesets = 'Cathode_top'
    new_block_name = 'cathode_s_block'
  [../]
  [./secondary_anode_block]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_boundary_cathode
    sidesets = 'Anode_bottom'
    new_block_name = 'anode_s_block'
  [../]
  [./primary_anode_block]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_anode_block
    sidesets = 'Electrolyte_top'
    new_block_name = 'anode_p_block'
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
    block = 'Anode Cathode Electrolyte'
    extra_vector_tags = 'ref'
    # use_finite_deform_jacobian = true
  [../]
[]

[Variables]
  [./ux]
    block = 'Anode Cathode Electrolyte'
  [../]
  [./uy]
    block = 'Anode Cathode Electrolyte'
  [../]
  [./V]
    block = 'Cathode Anode Electrolyte'
  [../]
  [./V_cathode]
    block = 'cathode_s_block'
  [../]
  [./V_anode]
    block = 'anode_s_block'
  [../]
  [./li_metal_conc]
    block = 'Cathode Anode'
  [../]
  [./normal_lm_x_cathode]
    block = 'cathode_s_block'
  [../]

  [./normal_lm_y_cathode]
    block = 'cathode_s_block'
  [../]

  [./normal_lm_x_anode]
    block = 'anode_s_block'
  [../]

  [./normal_lm_y_anode]
    block = 'anode_s_block'
  [../]

[]

[ICs]
  [./Voltage_cathode]
    type = ConstantIC
    block = 'Cathode'
    value = 3000.0
    variable = V
  [../]
  [./Voltage_electrolyte]
    type = ConstantIC
    block = 'Electrolyte Anode'
    value = 0.0
    variable = V
  [../]
  [./conc_cathode]
    type = ConstantIC
    block = 'Cathode'
    value = ${c_init_cat}
    variable = li_metal_conc
  [../]
  [./conc_anode]
    type = ConstantIC
    block = 'Anode'
    value = ${c_init_anode}
    variable = li_metal_conc
  [../]

[]

[Constraints]
  [./cathode_constraint]
    type = GapDisplacementConductanceConstraint
    variable = V_cathode
    secondary_variable = V
    secondary_boundary =  'Cathode_top'
    secondary_subdomain = 'cathode_s_block'
    primary_boundary = 'Electrolyte_bottom'
    primary_subdomain = 'cathode_p_block'
    # k_function = gapk
    k = 1e-3
    # use_displaced_mesh = true
    displacements = 'ux uy'
    compute_lm_residuals = true
    include_equilibrium_potential = true
    include_gap = false
    # extra_vector_tags = 'ref'
  [../]
  [./anode_constraint]
    type = GapDisplacementConductanceConstraint
    variable = V_anode
    secondary_variable = V
    secondary_boundary =  'Anode_bottom'
    secondary_subdomain = 'anode_s_block'
    primary_boundary = 'Electrolyte_top'
    primary_subdomain = 'anode_p_block'
    k_function = anode_rct
    k = ${k_anode}
    # use_displaced_mesh = true
    displacements = 'ux uy'
    compute_lm_residuals = true
    include_equilibrium_potential = false
    include_gap = true
    # extra_vector_tags = 'ref'
  [../]
  [./cathode_conc]
    type = ScaledBCConstraint
    variable = V_cathode
    secondary_variable = li_metal_conc
    primary_variable = V
    primary = false
    secondary_boundary =  'Cathode_top'
    secondary_subdomain = 'cathode_s_block'
    primary_boundary = 'Electrolyte_bottom'
    primary_subdomain = 'cathode_p_block'
    scale = -1.036426e-2
    extra_vector_tags = 'ref'
  [../]
  [./anode_conc]
    type = ScaledBCConstraint
    variable = V_anode
    secondary_variable = li_metal_conc
    primary_variable = V
    primary = false
    secondary_boundary =  'Anode_bottom'
    secondary_subdomain = 'anode_s_block'
    primary_boundary = 'Electrolyte_top'
    primary_subdomain = 'anode_p_block'
    scale = -1.036426e-2
    extra_vector_tags = 'ref'
  [../]
  [./cathode_mech_x]
    type = EqualValueConstraint
    variable = normal_lm_x_cathode
    primary_variable = ux
    secondary_variable = ux
    secondary_boundary =  'Cathode_top'
    secondary_subdomain = 'cathode_s_block'
    primary_boundary = 'Electrolyte_bottom'
    primary_subdomain = 'cathode_p_block'
    use_displaced_mesh = true
  [../]
  [./cathode_mech_y]
    type = EqualValueConstraint
    variable = normal_lm_y_cathode
    primary_variable = uy
    secondary_variable = uy
    secondary_boundary =  'Cathode_top'
    secondary_subdomain = 'cathode_s_block'
    primary_boundary = 'Electrolyte_bottom'
    primary_subdomain = 'cathode_p_block'
    use_displaced_mesh = true
  [../]
  [./anode_mech_x]
    type = EqualValueConstraint
    variable = normal_lm_x_anode
    primary_variable = ux
    secondary_variable = ux
    secondary_boundary =  'Anode_bottom'
    secondary_subdomain = 'anode_s_block'
    primary_boundary = 'Electrolyte_top'
    primary_subdomain = 'anode_p_block'
    use_displaced_mesh = true
  [../]
  [./anode_mech_y]
    type = EqualValueConstraint
    variable = normal_lm_y_anode
    primary_variable = uy
    secondary_variable = uy
    secondary_boundary =  'Anode_bottom'
    secondary_subdomain = 'anode_s_block'
    primary_boundary = 'Electrolyte_top'
    primary_subdomain = 'anode_p_block'
    use_displaced_mesh = true
  [../]


[]

[Functions]
  [./reaction_rate]
    type = PiecewiseLinear
    data_file = 'nmc_equilibrium_potential.csv'
    format = columns
  [../]
  [./current_density]
    type = ParsedFunction
    value = 'if (t < 10, 2e-3*t, 20e-3)'
  [../]
  [./anode_rct]
    type = ParsedFunction
    value = 'if (x < 3000, ${k_anode}, if (x > 4200, ${k_anode}, 0))'
  [../]
[]

[AuxVariables]
  [./flux_x]
    order = FIRST
    family = MONOMIAL
    block = 'Cathode Electrolyte Anode'
  [../]
  [./flux_y]
    order = FIRST
    family = MONOMIAL
    block = 'Cathode Electrolyte Anode'
  [../]
  [./flux_z]
    order = FIRST
    family = MONOMIAL
    block = 'Cathode Electrolyte Anode'
  [../]
  [./bndliflux]
    order = CONSTANT
    family = MONOMIAL
    block = 'Cathode Electrolyte Anode'
  [../]
  [./Eq_pot]
    order = CONSTANT
    family = MONOMIAL
    block = 'Cathode'
  [../]
  [./li_metal_flux_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./li_metal_flux_y]
    order = CONSTANT
    family = MONOMIAL
  [../]

[]

[AuxKernels]
  [./li_metal_flux_x]
    type = AnisoTropicDiffusionFluxAux
    diffusion_variable = li_metal_conc
    diffusivity = diffusivity_Li
    component = x
    block = 'Anode'
    variable = li_metal_flux_x
  [../]
  [./li_metal_flux_y]
    type = AnisoTropicDiffusionFluxAux
    diffusion_variable = li_metal_conc
    diffusivity = diffusivity_Li
    component = y
    block = 'Anode'
    variable = li_metal_flux_y
  [../]
  [./bnd_li_flux]
    type = DiffusionFluxNormalToBoundaryAux
    variable = bndliflux
    boundary = 'Electrolyte_bottom Cathode_top Electrolyte_top'
    diffusion_variable = V
    diffusivity = thermal_conductivity
  [../]
  [./Eq_pot]
    type = ADMaterialRealAux
    variable = Eq_pot
    property = equilibrium_potential
    block = 'Cathode'
  [../]


  [./li_ion_flux_x]
    type = ADDiffusionFluxAux
    variable = flux_x
    component = x
    diffusion_variable = V
    diffusivity = thermal_conductivity
    block = 'Cathode Electrolyte Anode'
  [../]

  [./li_ion_flux_y]
    type = ADDiffusionFluxAux
    variable = flux_y
    component = y
    diffusion_variable = V
    diffusivity = thermal_conductivity
    block = 'Cathode Electrolyte Anode'
  [../]
  [./li_ion_flux_z]
    type = ADDiffusionFluxAux
    variable = flux_z
    component = z
    diffusion_variable = V
    diffusivity = thermal_conductivity
    block = 'Cathode Electrolyte Anode'
  [../]
[]


[Kernels]
  [./li_ion_V]
    type = ADHeatConduction
    variable = V
    use_displaced_mesh = false
    block = 'Cathode Electrolyte Anode'
  [../]
  # [./li_metal2]
  #   type = ADChemoMechanoDiffusion
  #   variable = li_metal_conc
  #   diffusivity = diffusivity
  #   use_displaced_mesh = false
  #   block = 'Cathode'
  # [../]

  [./li_metal1]
    type = ADChemoMechanoAnsioDiffusion
    variable = li_metal_conc
    diffusivity = diffusivity_Li
    use_displaced_mesh = false
    block = 'Anode Cathode'
  [../]

  [./li_metal_dt]
    type = ADTimeDerivative
    variable = li_metal_conc
    use_displaced_mesh = false
    block = 'Cathode Anode'
  [../]

[]

[Materials]
  # [./diffusivity_Li]
  #   type = ADGenericConstantMaterial
  #   prop_names = 'diffusivity'
  #   prop_values = '5e8'
  #   block = 'Cathode'
  # [../]

  [./diffusivity_Li2]
    type = ADConstantAnisotropicMobility
    M_name = diffusivity_Li
    tensor = '0 0 0
              0 5e8 0
              0 0 0'
    block = 'Anode Cathode'
  [../]

  [./thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 100
    block = 'Cathode Anode'
  [../]
  [./thermal_conductivity2]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1
    block = 'Electrolyte'
  [../]

  [./equilibrium_potential]
    type = ADComputeEquilibriumPotential
    R = ${gas_constant}
    faraday = ${faraday}
    temperature = ${temperature}
    cref = ${c_ref_cat}
    concentration = li_metal_conc
    include_conc = false
    include_reaction_rate = true
    # reaction_rate = 780.0
    reaction_rate_function = reaction_rate
    include_mechanical_effects = false
    exclude_elastic_contribution = true
    block = 'Cathode'
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
    block = 'Electrolyte Cathode'
  [../]
  [./stress_llzo]
    type = ADComputeFiniteStrainElasticStress
    block = 'Cathode Electrolyte'
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
    gas_constant = ${gas_constant} # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = ${temperature}

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 2.0
    omega = 13
    alpha = '1 0 0'
    concentration = li_metal_conc
    cref = ${c_ref_anode}
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
    boundary = 'Cathode_bottom'
    # value = ${current_density}
    function = current_density
    variable = V
    extra_vector_tags = 'ref'
  [../]

  [./OV]
    type = ADDirichletBC
    boundary = 'Anode_top'
    value = 0.0
    variable = V
  [../]
  [./x_disp]
    type = ADDirichletBC
    value = 0
    boundary = 'Cathode_left Electrolyte_left Anode_left '
    variable = ux
    preset = true
  [../]
  [./y_disp]
    type = ADDirichletBC
    value = 0
    boundary = 'Cathode_bottom'
    variable = uy
    preset = true
  [../]

[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./Voltage_Cathode]
    type = SideAverageValue
    boundary = 'Cathode_bottom'
    variable = V
  [../]

  [./eq_pot]
    type = SideAverageValue
    boundary = 'Cathode_bottom'
    variable = Eq_pot
  [../]
  [./bot_conc]
    type = SideAverageValue
    boundary = 'Cathode_bottom'
    variable = li_metal_conc
  [../]
  [./top_conc]
    type = SideAverageValue
    boundary = 'Anode_top'
    variable = li_metal_conc
  [../]
  [./cathode_current]
    type = ADSideFluxIntegral
    variable = V
    diffusivity = thermal_conductivity
    boundary = 'Cathode_bottom'
  [../]
  [./anode_current]
    type = ADSideFluxIntegral
    variable = V
    diffusivity = thermal_conductivity
    boundary = 'Anode_top'
  [../]
  [./anode_metal_flux]
    type = SideAverageValue
    variable = li_metal_flux_y
    boundary = 'Anode_bottom'
  [../]
[]

[UserObjects]
  [./vcutoff]
    type = Terminator
    expression = 'Voltage_Cathode < ${v_cutoff}'
    execute_on = 'TIMESTEP_END'
  [../]
  [./max_conc]
    type = Terminator
    expression = 'bot_conc > 4.77e-2'
    execute_on = 'TIMESTEP_END'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  automatic_scaling = true
  # compute_scaling_once = false
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1     NONZERO               1e-15               '
  dt = 200
  # num_steps = 20
  nl_max_its = 25
  nl_abs_tol = 5e-10
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
  end_time = 50000
  # end_time =

  # num_steps =   3
  # snesmf_reuse_base = false
  compute_initial_residual_before_preset_bcs = true
  resid_vs_jac_scaling_param = 0.5
  verbose = true
  scaling_group_variables = 'V V_cathode V_anode'

[]
[Problem]
  type = ReferenceResidualProblem
  # solution_variables = 'ux uy normal_lm thermal_lm li_ion_V thermal_lm2 li_metal_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'V V_cathode V_anode li_metal_conc; ux uy normal_lm_x_cathode normal_lm_y_cathode normal_lm_x_anode normal_lm_y_anode'
  acceptable_iterations = 2
  coord_type = RZ
[]
[Outputs]
  # exodus = true
  # csv = true
  [./csv]
    type = CSV
    file_base = csv/test
  [../]
  [./out]
    type = Exodus
    file_base = rst/test
  [../]
[]
