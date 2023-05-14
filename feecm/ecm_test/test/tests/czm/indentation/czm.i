[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 5
    ny = 10
  []
  [subdomain_id]
    type = ParsedSubdomainMeshGenerator
    input = msh
    combinatorial_geometry = 'y > 0.5'
    block_id = 1
  []
  [split]
    type = BreakMeshByBlockGenerator
    input = subdomain_id
    split_interface = true
  []
  [right_corner]
    type = BoundingBoxNodeSetGenerator
    input = split
    bottom_left = '0.95 0.95 0'
    top_right = '1.05 1.05 0'
    new_boundary = 'right_corner'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Functions]
  [stretch]
    type = PiecewiseLinear
    x = '0 1 10 '
    y = '0 -1e3 -1e2'
  []
  [vel2]
    type = ParsedFunction
    value = '0.01'
  []
[]

[ICs]
  # [ux]
  #   type = RandomIC
  #   variable = disp_x
  #   min = 0
  #   max = 1e-3
  # []
  # [uy]
  #   type = RandomIC
  #   variable = disp_y
  #   min = 0
  #   max = 1e-3
  # []
[]

[BCs]
  [fix_x]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = 'left'
    variable = disp_x
  []
  [fix_y]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = 'bottom'
    variable = disp_y
  []
  [stretch]
    type = FunctionDirichletBC
    boundary = top
    variable = disp_y
    use_displaced_mesh = false
    function = 'stretch'
  []

  # [stretch]
  #   type = PresetVelocity
  #   boundary = top
  #   variable = disp_y
  #   use_displaced_mesh = false
  #   function = '0.0001'
  # []

[]

[Modules/TensorMechanics/CohesiveZoneMaster]
  [czm_ik_012]
    strain = SMALL
    boundary = 'Block0_Block1'
    generate_output = 'traction_x traction_y traction_z normal_traction tangent_traction jump_x '
                      'jump_y jump_z normal_jump tangent_jump'
    # base_name = 'czm_b012'
  []
[]

[Materials]
  # cohesive materials
  [czm_elastic_incremental]
    type = SalehaniIrani3DCTractionViscosity
    # type = SalehaniIrani3DCTraction
    # type = GaoBower3DCTraction
    boundary = 'Block0_Block1'
    maximum_normal_traction = 0.5e3
    maximum_shear_traction = 1e3
    normal_gap_at_maximum_normal_traction = 0.1
    tangential_gap_at_maximum_shear_traction = 0.1
    normal_viscosity = 0.000549
    tangential_viscosity = 0.000549
    # base_name = 'czm_b012'
  []
  # bulk materials
  [stress_0]
    type = ADComputeFiniteStrainElasticStress
    block = 1
  []
  [elasticity_tensor_0]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 200e3
    poissons_ratio = 0.33
    block = 1
  []
  [stress_1]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    perform_finite_strain_rotations = true
    block = 0
  []
  [elasticity_tensor_1]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 200e3
    poissons_ratio = 0.33
    block = 0
  []
  [plas]
    type = ADIsoTropicHyperVisco
    absolute_tolerance = 1e-6
    block = 0
    # relative_tolerance = 1e-06
    hardening_exponent = 1.0
    saturation_resistance = 2.0e3
    initial_resistance = 1e3
    hardening_modulus = 100.0
    rate_exponent = 0.05
    reference_strain_rate = 0.1
    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 1.0
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
  []
[]

[Modules]
  [TensorMechanics]
    [Master]
      [all]
        strain = FINITE
        add_variables = true
        # use_finite_deform_jacobian = true
        use_automatic_differentiation = true
        generate_output = 'stress_xx stress_yy stress_zz stress_xy'
      []
    []
  []
[]

[Postprocessors]
  [U]
    type = SideAverageValue
    variable = disp_y
    boundary = top
  []
  [sigma]
    type = SideAverageValue
    variable = stress_yy
    boundary = top
  []
  [jump]
    type = InterfaceAverageVariableValuePostprocessor
    boundary = Block0_Block1
    variable = jump_y
    interface_value_type = average
  []
  [normalized_U]
    type = ParsedPostprocessor
    function = 'U/0.1'
    pp_names = 'U'
  []
  [normalized_sigma]
    type = ParsedPostprocessor
    function = 'sigma/500'
    pp_names = 'sigma'
  []
  [normalized_jump]
    type = ParsedPostprocessor
    function = '2.0*jump/0.1'
    pp_names = 'jump'
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  # Executioner
  type = Transient

  solve_type = 'NEWTON'
  line_search = none
  petsc_options_iname = '-pc_type '
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-6
  l_max_its = 20
  start_time = 0.0
  dt = 0.01
  dtmax = 0.1
  dtmin = 0.000001
  end_time = 1
  # num_steps = 5
  automatic_scaling = true
  resid_vs_jac_scaling_param = 0.5
  compute_scaling_once = false
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-4
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 30
    # timestep_limiting_postprocessor = matl_ts_min
  []
[]

[Outputs]
  file_base = 'czm_bower_viscosity_0.1'
  exodus = true
  csv = true
[]
