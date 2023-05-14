# Gao Bower paper indentation verification with czm1, MSMSE, 2004
# Geometry and bcs from Abdul Baqi and Van der Giessen JMR 2001

[Mesh]
  patch_update_strategy = iteration
  patch_size = 200
  parallel_type = REPLICATED
  [mesh]
    type = FileMeshGenerator
    file = bimat2.e
  []
  [split]
    type = BreakMeshByBlockGenerator
    input = mesh
    split_interface = true
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Functions]
  [indenter]
    type = PiecewiseLinear
    x = '0 1 10'
    y = '0 -1e-1 -1'
  []
[]

[Modules]
  [TensorMechanics]
    [Master]
      [all]
        strain = FINITE
        add_variables = true
        use_finite_deform_jacobian = true
        use_automatic_differentiation = true
        generate_output = 'stress_xx stress_yy stress_xy vonmises_stress strain_xx strain_yy '
                          'strain_xy'
        # block = 'Film Substrate'
      []
    []
  []
[]

[Modules/TensorMechanics/CohesiveZoneMaster]
  [czm]
    strain = SMALL
    boundary = 'Substrate_Film'
    generate_output = 'traction_x traction_y traction_z normal_traction tangent_traction jump_x '
                      'jump_y jump_z normal_jump tangent_jump'
  []
[]

[Materials]
  # cohesive materials
  [czm_elastic_incremental]
    # type = SalehaniIrani3DCTractionViscosity
    type = GaoBower3DCTraction
    boundary = 'Substrate_Film'
    maximum_normal_traction = 500
    maximum_shear_traction = 1000
    normal_gap_at_maximum_normal_traction = 0.1
    tangential_gap_at_maximum_shear_traction = 0.1
    normal_viscosity = 52.5e-6
    tangential_viscosity = 52.5e-6
    # base_name = 'czm_b012'
  []

  [stress_film_indenter]
    type = ADComputeLinearElasticStress
    block = 'Film'
  []
  [elasticity_tensor_film]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 500e3
    poissons_ratio = 0.33
    block = 'Film'
  []

  [elasticity_tensor_substrate]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 200e3
    poissons_ratio = 0.33
    block = 'Substrate'
  []
  [stress_substrate]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    block = 'Substrate'
    perform_finite_strain_rotations = true
  []
  [plas]
    type = ADIsoTropicHyperVisco
    block = 'Substrate'
    absolute_tolerance = 1e-6
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

[Problem]
  coord_type = RZ
[]

[BCs]
  [axisymmetry]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = 'block_left'
  []
  [bottom]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = 'block_bottom'
  []
  [indenter]
    type = FunctionDirichletBC
    variable = disp_y
    function = 'indenter'
    boundary = 'Film_top'
  []
[]

[Outputs]
  exodus = true
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
  dt = 1e-3
  dtmax = 0.1
  dtmin = 1e-6
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
