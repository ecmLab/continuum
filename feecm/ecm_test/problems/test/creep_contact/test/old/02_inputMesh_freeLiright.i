order = FIRST
strnRt = 0.06
name = '02_inputMesh_freeLiright'

[Mesh]
  patch_size = 80
  patch_update_strategy = auto
  [./blockSE]
    type = FileMeshGenerator
    file = data/blockSE.msh
  [../]
  [./blockSE_sidesets]
    type = SideSetsFromNormalsGenerator
    input = blockSE
    normals = '0  -1  0
               1   0  0
               0   1  0
              -1   0  0'
    fixed_normal = true
    new_boundary = 'blockSE_bottom blockSE_right blockSE_top blockSE_left'
  [../]
  [./blockSE_id]
    type = SubdomainIDGenerator
    input = blockSE_sidesets
    subdomain_id = 1
  [../]

  [./blockLi]
    type = FileMeshGenerator
    file = data/blockLi.msh
  [../]
  [./blockLi_id]
    type = SubdomainIDGenerator
    input = blockLi
    subdomain_id = 2
  [../]

  [./combined]
    type = MeshCollectionGenerator
    inputs = 'blockSE_id blockLi_id'
  [../]
  [./block_rename]
    type = RenameBlockGenerator
    input = combined
    old_block_id = '1 2'
    new_block_name = 'blockSE blockLi'
  [../]
  [./block_sidesets]
    type = SideSetsFromPointsGenerator
    input = block_rename
    points = '0.0   0.0001   0
              0.25   0.25  0
              0.0   0.5   0
              0.0  0.25  0'
    new_boundary = 'blockLi_bottom blockLi_right blockLi_top blockLi_left'
  [../]

[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'disp_x disp_y'
  acceptable_iterations = 2
[]

[Variables]
  [./disp_x]
   order = ${order}
   block = 'blockSE blockLi'
  [../]
  [./disp_y]
   order = ${order}
   block = 'blockSE blockLi'
  [../]
[]

[Functions]
  [./pressure]
    type = ParsedFunction
    value = 'if (t <= 5, 0.24*t, 1.2)'
  [../]
  [./uy]
    type = ParsedFunction
    value = 'if (t <=1, 0.05*t, 0.05)'
  [../]
  [./vel2]
    type = ParsedFunction
    value  = 'if (t <= 1, ${strnRt} * exp(${strnRt}*t),0)'
  [../]
[]

[Modules/TensorMechanics/Master]
  [./action]
    add_variables = true
    strain = FINITE
    use_displaced_mesh = true
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy vonmises_stress strain_xx strain_yy'
    block = 'blockSE blockLi'
    extra_vector_tags = 'ref'
    use_automatic_differentiation = true
  [../]
[]

[Contact]
  [./frictionless]
#    mesh = block_rename
    mesh = block_sidesets
    primary = blockSE_top
    secondary = blockLi_bottom
    penalty = 1e7
#    formulation = kinematic
    formulation = mortar
  [../]
[]

[BCs]
#  [./pressure]
#    type = ADPressure
#    variable = disp_y
#    component = 1
#    function = pressure
#    boundary = 'blockSE_bottom'
#  [../]
#  [./uy]
#    type = ADFunctionDirichletBC
#    variable = disp_y
#    function = uy
#    boundary = 'blockSE_bottom'
#  [../]
  [./vely]
    type = PresetVelocity
    variable = disp_y
    boundary = 'blockSE_bottom'
    function = vel2
  [../]
  [./left_right]
    type = ADDirichletBC
    variable = disp_x
    boundary = 'blockLi_left blockSE_left blockSE_right blockSE_bottom'
    value = 0
    preset = true
  [../]
  [./bottom]
    type = ADDirichletBC
    variable = disp_y
    boundary = 'blockLi_top'
    value = 0
    preset = true
  [../]

[]

[Materials]
  [./blockSE]
    type = ADComputeIsotropicElasticityTensor
    block = 'blockSE'
    poissons_ratio = 0.3
    youngs_modulus = 130e9
  [../]
  [./blockLi]
    type = ADComputeIsotropicElasticityTensor
    block = 'blockLi'
    poissons_ratio = 0.38
    youngs_modulus = 7.8e9
  [../]
  [./stress_blockSE]
    type = ADComputeFiniteStrainElasticStress
    block = 'blockSE'
  [../]
 
  [./stress_blockLi]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    # perform_finite_strain_rotations = true
    block = 'blockLi'
  [../]
  [./plas]
    type = ADIsoTropicHyperViscoCreep
    # absolute_tolerance = 1e-6
    # relative_tolerance = 1e-06
    hardening_exponent = 2.0
    saturation_resistance = 2.0e6
    initial_resistance = 0.95e6
    hardening_modulus = 10.0e6
    rate_exponent = 0.15
    activation_energy = 37
    gas_constant = 8.314462681e-3 # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = 298

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 1.0
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    block = 'blockLi'
  [../]

[]

[Postprocessors]
  [./stress_yy]
    type = SideAverageValue
    variable = stress_yy
    boundary = blockLi_top
  [../]
  [./matl_ts_min]
    type = MaterialTimeStepPostprocessor
    block = 'blockLi'
  [../]
 [./contact]
    type = ContactDOFSetSize
    variable = frictionless_normal_lm
    subdomain = frictionless_secondary_subdomain
  [../]
[]

[Executioner]
  type = Transient
## Solver setting
  solve_type = 'NEWTON'
  automatic_scaling = true
  petsc_options = '-snes_converged_reason -ksp_converged_reason'
  petsc_options_iname = '-pc_type -mat_mffd_err -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu       1e-5          NONZERO               1e-15'

#  solve_type = 'PJFNK'
#  automatic_scaling = true
#  petsc_options = '-snes_ksp_ew'
#  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
#  petsc_options_value = 'hypre    boomeramg      101'

  line_search = 'contact'
  nl_abs_tol = 1e-10
  l_max_its = 100
  nl_max_its = 20

## Timestep setting
  dt = 0.005
  dtmin = 1e-5
  dtmax = 0.5
  end_time = 20
  timestep_tolerance = 1e-6
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.005
    growth_factor = 1.5
    cutback_factor = 0.65
    optimal_iterations = 40
    timestep_limiting_postprocessor = matl_ts_min
  [../]

[]

[Outputs]
  exodus = true
  file_base = rst/${name}
  [out]
    type = Checkpoint
    interval = 5
    num_files = 2
  []
[]

[Debug]
  show_var_residual_norms = true
[]
