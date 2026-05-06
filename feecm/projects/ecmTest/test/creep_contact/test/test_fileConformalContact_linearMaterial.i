order = FIRST
strnRt = 0.02
name = 'test_linearMaterial'

[Mesh]
  patch_size = 80
  patch_update_strategy = auto
  [./blocks]
    type = FileMeshGenerator
    file = data/conformalContact.msh
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
   block = 'blockCeramic blockMetal'
  [../]
  [./disp_y]
   order = ${order}
   block = 'blockCeramic blockMetal'
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
    block = 'blockCeramic blockMetal'
    extra_vector_tags = 'ref'
    use_automatic_differentiation = true
  [../]
[]

[Contact]
  [./frictionless]
    mesh = blocks
    primary = blockCeramic_top
    secondary = blockMetal_bottom
    penalty = 1e5
#    formulation = kinematic
    formulation = mortar
  [../]
[]

[BCs]
  [./vely]
    type = PresetVelocity
    variable = disp_y
    boundary = 'blockCeramic_bottom'
    function = vel2
  [../]
  [./left_right]
    type = ADDirichletBC
    variable = disp_x
    boundary = 'blockMetal_left blockCeramic_left blockCeramic_right blockCeramic_bottom'
    value = 0
    preset = true
  [../]
  [./bottom]
    type = ADDirichletBC
    variable = disp_y
    boundary = 'blockMetal_top'
    value = 0
    preset = true
  [../]

[]

[Materials]
  [./blockCeramic]
    type = ADComputeIsotropicElasticityTensor
    block = 'blockCeramic'
    poissons_ratio = 0.3
    youngs_modulus = 130e9
  [../]
  [./blockMetal]
    type = ADComputeIsotropicElasticityTensor
    block = 'blockMetal'
    poissons_ratio = 0.38
    youngs_modulus = 7.8e9
  [../]
  [./stress_blockCeramic]
    type = ADComputeFiniteStrainElasticStress
    block = 'blockCeramic blockMetal'
  [../]

[]

[Postprocessors]
  [./stress_yy]
    type = SideAverageValue
    variable = stress_yy
    boundary = blockMetal_top
  [../]
#  [./matl_ts_min]
#    type = MaterialTimeStepPostprocessor
#    block = 'blockMetal'
#  [../]
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
  end_time = 5
  timestep_tolerance = 1e-6
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.005
    growth_factor = 1.5
    cutback_factor = 0.65
    optimal_iterations = 40
#    timestep_limiting_postprocessor = matl_ts_min
  [../]

[]

[Outputs]
  exodus = true
  file_base = rst/${name}
[]

[Debug]
  show_var_residual_norms = true
[]
