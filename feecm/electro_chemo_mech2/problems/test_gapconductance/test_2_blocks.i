[Mesh]
  patch_size = 80
  patch_update_strategy = auto
  parallel_type = REPLICATED
  # file = check/full_model_cp/0002_mesh.cpr
  [./mesh]
    type = FileMeshGenerator
    file = data/2blocks.e
  [../]
  [./primary_block]
    type = LowerDBlockFromSidesetGenerator
    input = mesh
    sidesets = 'top_bottom'
    new_block_name = 'mech_contact_secondary_subdomain'
  [../]
  [./secondary_block]
    type = LowerDBlockFromSidesetGenerator
    input = primary_block
    sidesets = 'bottom_top'
    new_block_name = 'mech_contact_primary_subdomain'
  [../]

[]

[Functions]
  [./pressure]
    type = ParsedFunction
    # value = 'if (t<=1,0.2 + 1.8*t,2.0)'
    value = 'if (t<=1,0.01, 0.01)'
    # type = PiecewiseLinear
    # data_file = 'pressure_time1.csv'
    # format = columns
  [../]
  [./flux]
    type = ParsedFunction
    value = 'if (t < 1.0, 1e-3*(t), 1e-3)'
    # value = '1e-4'
  [../]
  [./gapk]
    type = PiecewiseLinear
    x = '0 1e-15 1e-14 1    1e6   1e9'
    y = '0  0    1e-3  1e-2 1e-3    1e-3 '
  [../]
  [./gapk1]
    type = PiecewiseLinear
    x = '-5   -1e-3 -1e-6 -1e-15 1e-12 1e-11 1'
    y = '1e-3  1e-3  1e-3   1e-3  1e-3  0    0'
  [../]

[]

[GlobalParams]
  displacements = 'ux uy'
[]

[Variables]
  [./ux]
    block = 'blockTop blockBottom'
  [../]
  [./uy]
    block = 'blockTop blockBottom'
  [../]
  [./li_ion_V]
    block = 'blockTop blockBottom'
  [../]
  [./thermal_lm]
    block = 'mech_contact_secondary_subdomain'
  [../]
  [./normal_lm_x]
    block = 'mech_contact_secondary_subdomain'
  [../]

  [./normal_lm_y]
    block = 'mech_contact_secondary_subdomain'
  [../]

  [./li_metal_conc]
    block = 'blockTop'
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress hydrostatic_stress'
    use_automatic_differentiation = true
    extra_vector_tags = 'ref'
    # use_finite_deform_jacobian = true
    block = 'blockTop blockBottom'
  [../]
[]


[Problem]
  type = ReferenceResidualProblem
  # solution_variables = 'ux uy normal_lm thermal_lm li_ion_V thermal_lm2 li_metal_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'ux uy; li_ion_V li_metal_conc thermal_lm'
  acceptable_iterations = 2
  # restart_file_base = check/full_model_cp/0002
[]

[Constraints]
  [./thermal_constraint]
    type = GapDisplacementConductanceConstraint
    # contact_pressure = mech_contact_normal_lm
    variable = thermal_lm
    secondary_variable = li_ion_V
    primary_boundary =  'bottom_top'
    secondary_subdomain = 'mech_contact_secondary_subdomain'
    secondary_boundary = 'top_bottom'
    primary_subdomain = 'mech_contact_primary_subdomain'
    k = 1e-3
    use_displaced_mesh = true
    include_gap = false
    include_equilibrium_potential = false
    R = 8.3145
    faraday = 96.4853329
    temperature = 298
    surfaceType = SECONDARY
    include_concentration = true
    concentration = li_metal_conc
    # compute_lm_residuals = true
    # compute_primal_residuals = true
  [../]

  [./therm2]
    type = ScaledBCConstraint
    variable = thermal_lm
    primary_variable = li_ion_V
    secondary_variable = li_metal_conc
    primary_boundary =  'bottom_top'
    secondary_subdomain = 'mech_contact_secondary_subdomain'
    secondary_boundary = 'top_bottom'
    primary_subdomain = 'mech_contact_primary_subdomain'
    scale = 1.036426e-2
    use_displaced_mesh = true
    primary = false
  [../]

    [./normal_lm_x]
      type = EqualValueConstraint
      variable = normal_lm_x
      primary_variable = ux
      secondary_variable = ux
      primary_boundary =  'bottom_top'
      secondary_subdomain = 'mech_contact_secondary_subdomain'
      secondary_boundary = 'top_bottom'
      primary_subdomain = 'mech_contact_primary_subdomain'
      # delta = 1.0
      # delta = 0.1
    [../]

    [./normal_lm_y]
      type = EqualValueConstraint
      variable = normal_lm_y
      primary_variable = uy
      secondary_variable = uy
      primary_boundary =  'bottom_top'
      secondary_subdomain = 'mech_contact_secondary_subdomain'
      secondary_boundary = 'top_bottom'
      primary_subdomain = 'mech_contact_primary_subdomain'
      # delta = 0.5
    [../]

[]


[Kernels]
  [./li_ion_V]
    type = ADHeatConduction
    block = 'blockTop blockBottom'
    variable = li_ion_V
    use_displaced_mesh = false
  [../]
  [./li_metal2]
    type = ADMatAnisoDiffusion
    block = 'blockTop'
    variable = li_metal_conc
    diffusivity = diffusivity
    use_displaced_mesh = false
  [../]
  [./li_metal_dt]
    type = ADTimeDerivative
    block = 'blockTop'
    variable = li_metal_conc
    use_displaced_mesh = false
  [../]

[]


[Materials]
  [./diffusivity_Li]
    type = ADDiffusionAlongPrincipalDirectionsMaterial
    diffusivity_vector = '1e8 0 0'
    block = 'blockTop'
  [../]
  [./thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1.0
    block = 'blockBottom'
  [../]

  [./thermal_conductivity2]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1e2
    block = 'blockTop'
    use_displaced_mesh = true
  [../]

  [./elasticity_tensor_Li]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 7.8e3
    poissons_ratio = 0.3
    block = 'blockTop'
  [../]
  [./elasticity_tensor_Llzo]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 170e3
    poissons_ratio = 0.25
    block = 'blockBottom'
  [../]
  [./stress_llzo]
    type = ADComputeFiniteStrainElasticStress
    block = 'blockBottom'
  [../]

  [./stress_Li]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    # perform_finite_strain_rotations = true
    block = 'blockTop'
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
    block = 'blockTop'
    intBnd = 'bottom_top'
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    enable = true
    output_properties = 'growth_stress growth_tensor'
  [../]

[]

[BCs]
  [./Li_top_y]
    type = ADDirichletBC
    variable = uy
    boundary = 'bottom_bottom'
    value = 0.0
  [../]


  [./left_right]
    type = ADDirichletBC
    variable = ux
    boundary = 'bottom_left top_left'
    value = 0
  [../]
  [./li_ion_bottom]
    type = FunctionNeumannBC
    boundary = 'bottom_bottom'
    variable = li_ion_V
    function = flux
    extra_vector_tags = 'ref'
  [../]

  [./li_ion_top]
    type = ADDirichletBC
    boundary = 'top_top'
    variable = li_ion_V
    value = 0.0
  [../]

  [./top_pressure]
    type = ADPressure
    boundary = 'top_top'
    variable = uy
    component = 1
    function = pressure
    use_displaced_mesh = true
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
  # compute_scaling_once = false
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor'
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_force_iteration -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist 1     NONZERO               1e-20               '
  # line_search = basic
  line_search = contact
  # line_search = none
  l_max_its = 50
  nl_max_its = 15
  nl_abs_tol = 1e-14
  nl_rel_tol = 1e-3
  # start_time = 0.0
  dt = 1.0
  dtmax = 5.0
  dtmin = 1e-5
  # end_time = 5000.0
    num_steps = 5
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 40
    # timestep_limiting_postprocessor = matl_ts_min
  [../]
  resid_vs_jac_scaling_param = 0.5
  scaling_group_variables = 'li_ion_V thermal_lm'
[] # Executioner

[Outputs]
  [./out]
    type = Exodus
    interval = 2
    file_base = rst/full_model
    # sync_times = 1.0
  [../]
  [./csv]
    type = CSV
    file_base = csv/full_model
  [../]
   # [./checkpoint]
   #   type = Checkpoint
   #   file_base = check/full_model
   #   num_files = 4
   #   interval = 20
   # [../]

[] # Outputs

[Debug]
  show_var_residual_norms = true
  # show_material_props = true
[]

[Postprocessors]
  [./bottom_current]
    type = ADSideFluxIntegral
    variable = li_ion_V
    boundary = 'bottom_top'
    diffusivity = thermal_conductivity
  [../]
  [./top_current]
    type = ADSideFluxIntegral
    variable = li_ion_V
    boundary = 'top_bottom'
    diffusivity = thermal_conductivity
  [../]
  [./lagrange]
    type = ElementIntegralVariablePostprocessor
    variable = thermal_lm
    block = 'mech_contact_secondary_subdomain'
  [../]
  # [./auxcurrent]
  #   type = SideIntegralVariablePostprocessor
  #   variable = bndliflux
  #   boundary = 15
  # [../]
  [./over_potential]
    type = SideAverageValue
    variable = li_ion_V
    boundary = 'bottom_top'
  [../]

  [./ext_pressure]
    type = SideAverageValue
    variable = stress_yy
    boundary = 'bottom_bottom'
  [../]
  [./matl_ts_min]
    type = MaterialTimeStepPostprocessor
    block = 'blockTop'
  [../]
  # [./contact_size]
  #   type = ContactDOFSetSize
  #   variable = mech_contact_normal_lm
  #   subdomain = 'mech_contact_secondary_subdomain'
  #   tolerance = 1e-15
  # [../]

  [./limetal_conc]
    type = SideAverageValue
    variable = li_metal_conc
    boundary = 'top_bottom'
  [../]
[]
