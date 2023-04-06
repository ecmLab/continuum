[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 1
    ny = 2
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
    x = '0 1'
    y = '0 0.1'
  []
[]

[BCs]
  [fix_x]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = bottom
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
    function = stretch
  []

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
    boundary = 'Block0_Block1'
    maximum_normal_traction = 182
    maximum_shear_traction = 364
    normal_gap_at_maximum_normal_traction = 0.01
    tangential_gap_at_maximum_shear_traction = 0.01
    normal_viscosity = 0.0549
    tangential_viscosity = 0.054
    # base_name = 'czm_b012'
  []
  # bulk materials
  [stress]
    type = ADComputeFiniteStrainElasticStress
  []
  [elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 200e4
    poissons_ratio = 0.3
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
    function = 'U/0.01'
    pp_names = 'U'
  []
  [normalized_sigma]
    type = ParsedPostprocessor
    function = 'sigma/182'
    pp_names = 'sigma'
  []
  [normalized_jump]
    type = ParsedPostprocessor
    function = '2.0*jump/0.01'
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
  # automatic_scaling = true
  # resid_vs_jac_scaling_param = 0.5
  # compute_scaling_once = false
[]

[Outputs]
  file_base = 'czm_bower_viscosity_0.1'
  exodus = true
  csv = true
[]
