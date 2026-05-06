
current_density = 10e-3
k_anode = 1e-3
faraday = 96.4853329
temperature = 298
gas_constant = 8.31446
c_ref = 0
c_ref_pore = 2.5925e-2
k_pore = 1e-3
c_init_Li = 1e-10
c_init_C = 1e-10
c_max_C = 3.05e-2
pressure = 1e-4
end_time = 150
[Mesh]
  [./mesh]
    type = FileMeshGenerator
    file = 'data/test_2_blocks3.e'
  [../]
[]
[Variables]
  [./ux]
    block = 'left1 right1'
  [../]
  [./uy]
    block = 'left1 right1'
  [../]
  [./li_conc]
    block = 'left1 right1'
  [../]
  [./equil]
    block = 'left_right_secondary_subdomain'
  [../]

[]

[ICs]
  [./li_conc_left]
    type = ConstantIC
    value = ${c_init_Li}
    variable = li_conc
    block = 'left1'
  [../]
  [./li_conc_right]
    type = ConstantIC
    value = ${c_init_C}
    variable = li_conc
    block = 'right1'
  [../]
[]

[GlobalParams]
  displacements = 'ux uy'
[]

[Contact]
  [./left_right]
      formulation = mortar
      normal_smoothing_distance = 0.5
      tangential_tolerance = 0.25
      primary = 'right_left'
      secondary = 'left_right'
      mesh = mesh
  [../]
[]
[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress hydrostatic_stress'
    use_automatic_differentiation = true
    block = 'left1 right1'
    extra_vector_tags = 'ref'
  [../]
[]
[AuxVariables]
  [./li_metal_flux_x]
    order = CONSTANT
    family = MONOMIAL
    block = 'right1'
  [../]
  [./li_metal_flux_y]
    order = CONSTANT
    family = MONOMIAL
    block = 'right1'
  [../]
  [./li_metal_flux_z]
    order = CONSTANT
    family = MONOMIAL
    block = 'right1'
  [../]
  [./eq_pot]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./bndliflux]
    order = CONSTANT
    family = MONOMIAL
    block = 'left1 right1'
  [../]
  [./normal_stress]
    order = FIRST
    family = MONOMIAL
    block = 'interLa'
  [../]
[]
[AuxKernels]
  [./bnd_li_flux]
    type = DiffusionFluxNormalToBoundaryAux
    variable = bndliflux
    boundary = 'right_left'
    diffusion_variable = li_conc
    diffusivity = diffusivity
  [../]
  [./li_metal_flux_x1]
    type = ADDiffusionFluxAux
    variable = li_metal_flux_x
    component = x
    diffusion_variable = li_conc
    block = 'right1'
    diffusivity = diffusivity
  [../]

  [./li_metal_flux_y1]
    type = ADDiffusionFluxAux
    variable = li_metal_flux_y
    component = y
    diffusion_variable = li_conc
    block = 'right1'
    diffusivity = diffusivity
  [../]
  [./li_metal_flux_z1]
    type = ADDiffusionFluxAux
    variable = li_metal_flux_z
    component = z
    diffusion_variable = li_conc
    block = 'right1'
    diffusivity = diffusivity
  [../]
  [./eq_pot]
    type = ADMaterialRealAux
    variable = eq_pot
    block = 'left1 right1'
    property = 'chemical_potential'
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
  [./equil]
    type = ParsedFunction
    value = 'if (t < 1e-5, 1e-16, 0)'
  [../]
[]



[Constraints]
  # [./left_right_conc]
  #   type = GapEquilibriateConstraint
  #   include_gap = false
  #   displacements = 'ux uy'
  #   k = -1e-3
  #   k_function = 'equil'
  #   variable = 'equil'
  #   use_displaced_mesh = true
  #   secondary_mat_prop = chemical_potential
  #   primary_mat_prop = chemical_potential
  #   primary_boundary = 'right_left'
  #   secondary_boundary = 'left_right'
  #   primary_subdomain = 'left_right_primary_subdomain'
  #   secondary_subdomain = 'left_right_secondary_subdomain'
  #   primary_variable = 'li_conc'
  #   secondary_variable = 'li_conc'
  #   extra_vector_tags = 'ref'
  #   faraday = ${faraday}
  #   temperature = ${temperature}
  #   R = ${gas_constant}
  #   prefactor = RT
  #   one_sided = NONE
  #   # one_sided = SECONDARY->PRIMARY
  # [../]
  [./left_right_conc]
    type = EqualValueConstraint
    primary_boundary = 'right_left'
    secondary_boundary = 'left_right'
    primary_subdomain = 'left_right_primary_subdomain'
    secondary_subdomain = 'left_right_secondary_subdomain'
    component = 0
    primary_variable = eq_pot
    secondary_variable = eq_pot
    variable = equil
  [../]
[]

[Kernels]
  [./li_conc]
    type = ADChemoMechanoDiffusion
    variable = li_conc
    diffusivity = diffusivity
    use_displaced_mesh = false
    block ='left1 right1'
  [../]

  [./li_metal_dt]
    type = ADTimeDerivative
    variable = li_conc
    use_displaced_mesh = false
    block = 'left1 right1'
  [../]
[]

[Materials]
  [./diffusivity]
    type = ADIsotropicDiffusionMaterial
    diffusion_coef = 1e6
    block = 'left1'
  [../]
  [./diff]
    type = ADIsotropicDiffusionMaterial
    diffusion_coef = 1.0
    block = 'right1'
  [../]
  [./equilibrium_potential_Li]
    type = ADComputeEquilibriumPotential
    R = ${gas_constant}
    faraday = ${faraday}
    temperature = ${temperature}
    cref = 1.0
    concentration = li_conc
    include_conc = true
    include_reaction_rate = true
    reaction_rate = 0.0
    is_ocv = true
    # reaction_rate_function = reaction_rate
    include_mechanical_effects = false
    exclude_elastic_contribution = true
    block = 'left1'
  [../]

  [./equilibrium_potential_pore]
    type = ADComputeEquilibriumPotential
    R = ${gas_constant}
    faraday = ${faraday}
    temperature = ${temperature}
    cref = 1.0
    concentration = li_conc
    include_conc = true
    include_reaction_rate = true
    reaction_rate = 0
    # reaction_rate_function = reaction_rate_c
    include_mechanical_effects = false
    exclude_elastic_contribution = true
    block = 'right1'
    potential = BOWER
  [../]

  [./elasticity_tensor_C]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 20e3
    poissons_ratio = 0.3
    block = 'right1'
  [../]

  [./Stress_tensor_C]
    type = ADComputeFiniteStrainElasticStress
    block = 'right1'
  [../]

  [./elasticity_tensor_Li]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 7.8e3
    poissons_ratio = 0.3
    block = 'left1'
  [../]
  [./stress_Li]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    # perform_finite_strain_rotations = true
    block = 'left1'
  [../]
  # [./stress_Li2]
  #   type = ADComputeMultipleInelasticStress
  #   inelastic_models = 'elastic'
  #   # perform_finite_strain_rotations = true
  #   block = 'right1'
  # [../]
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
    omega = 0.0
    alpha = '1 0 0'
    concentration = li_conc
    cref = ${c_init_Li}
    block = 'left1'
    intBnd = 'left_top'
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    enable = true
    output_properties = 'growth_stress'
    outputs = out
  [../]
  # [./elastic]
  #   type = ADIsotropicElasticSwelling
  #   # rate_form = false
  #   isotropic_swelling = true
  #   absolute_tolerance = 1e-5
  #   relative_tolerance = 1e-6
  #   alpha = '0.333333 0.333333 0.333333'
  #   # reference_strain_rate = 0.05
  #   temperature = ${temperature}
  #   omega = 0.1
  #   block = 'right1'
  #   internal_solve_full_iteration_history = true
  #   internal_solve_output_on = on_error
  #   enable = true
  #   cref = ${c_init_C}
  #   concentration = li_conc
  # [../]
[]

[BCs]
  [./ux]
    type = ADDirichletBC
    preset = true
    boundary = 'left_left right_right'
    variable = 'ux'
    value = 0
  [../]
  [./uy]
    type = ADDirichletBC
    preset = true
    boundary = 'left_bottom right_bottom'
    variable = 'ux'
    value = 0
  [../]
  # [./pressure]
  #   type = ADPressure
  #   variable = uy
  #   component = 1
  #   boundary = 'left_top right_top'
  #   function = pressure
  # [../]
  [./conc_li_flux]
    type = ADNeumannBC
    # boundary = 'left_left'
    boundary = 'right_right'
    variable = li_conc
    value = '1e-3'
    extra_vector_tags = 'ref'
  [../]
  # [./conc_li]
  #   type = ADDirichletBC
  #   value = ${c_init_Li}
  #   boundary = 'left_top left_left'
  #   variable = li_conc
  # [../]
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
  nl_abs_tol = 1e-9
  nl_rel_tol = 1e-3
  dtmax = 200
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-2
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 15
    # timestep_limiting_postprocessor = matl_ts_min
  [../]
  # num_steps = 2
  end_time = ${end_time}
  # scaling_group_variables = 'V flux_interlayer flux_pore'
  resid_vs_jac_scaling_param = 0.5
[]

[Problem]
  type = ReferenceResidualProblem
  # solution_variables = 'ux uy normal_lm thermal_lm V thermal_lm2 li_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'ux uy left_right_normal_lm; li_conc equil'
  acceptable_iterations = 2
  # coord_type = RZ
[]

[Postprocessors]
  [./equilibrating_flux]
    type = ElementAverageValue
    block = 'left_right_secondary_subdomain'
    variable = equil
  [../]
  [./Li_conc_li]
    type = SideAverageValue
    boundary = 'left_right'
    variable = li_conc
  [../]
  [./Li_conc_C]
    type = SideAverageValue
    boundary = 'right_left'
    variable = li_conc
  [../]
  [./chem_pot_C]
    type = SideAverageValue
    boundary = 'right_left'
    variable = eq_pot
  [../]

[]

[Outputs]
  # exodus = true
  # csv = true
  [./csv]
    type = CSV
    file_base = csv/test_pore_equil
  [../]
  [./out]
    type = Exodus
    file_base = rst/test_pore_equil
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
[]
