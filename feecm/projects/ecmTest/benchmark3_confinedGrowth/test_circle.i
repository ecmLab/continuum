order = FIRST
strnRt = 0.002
name = 'circle'

[Mesh]
  patch_size = 80
  patch_update_strategy = auto
  [./blocks]
    type = FileMeshGenerator
    file = data/circle_conformal.msh
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
   block = 'blockCeramic blockMetal interLayer'
  [../]
  [./disp_y]
   order = ${order}
   block = 'blockCeramic blockMetal interLayer'
  [../]
[]

[AuxVariables]
  [swellNormal]
    family = MONOMIAL_VEC
    block = interLayer
  []

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
    block = 'blockCeramic blockMetal interLayer'
    extra_vector_tags = 'ref'
    use_automatic_differentiation = true
  [../]
[]

[Contact]
  [./frictionless]
    mesh = blocks
    primary = blockCeramic_top
    secondary = blockMetal_bottom
    penalty = 1e-2
#    formulation = kinematic
    formulation = mortar
  [../]
[]

[AuxKernels]
  [normals]
    type = ADVectorMaterialRealVectorValueAux
    property = swell_normal
    variable = swellNormal
    block = interLayer
  []

[]

[BCs]
#  [./pressure]
#    type = ADPressure
#    variable = disp_y
#    component = 1
#    function = pressure
#    boundary = 'blockCeramic_bottom'
#  [../]
#  [./uy]
#    type = ADFunctionDirichletBC
#    variable = disp_y
#    function = uy
#    boundary = 'blockCeramic_bottom'
#  [../]
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
    youngs_modulus = 130e-3
  [../]
  [./blockMetal]
    type = ADComputeIsotropicElasticityTensor
    block = 'blockMetal interLayer'
    poissons_ratio = 0.38
    youngs_modulus = 7.8e-3
  [../]

  [./stress_blockCeramic]
    type = ADComputeFiniteStrainElasticStress
    block = 'blockCeramic'
  [../]
  [./stress_blockMetal]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'creep'
    # perform_finite_strain_rotations = true
    block = 'blockMetal'
  [../]
  [./stress_interLayer]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'creepSwell'
    # perform_finite_strain_rotations = true
    block = 'interLayer'
  [../]

  [./creep]
    type = ADIsoTropicHyperViscoCreep
    # absolute_tolerance = 1e-6
    # relative_tolerance = 1e-06
    hardening_exponent = 2.0
    rate_exponent = 0.15
    saturation_resistance = 2.0e-6
    hardening_modulus = 10.0e-6
    initial_resistance = 0.95e-6
    activation_energy = 37
    gas_constant = 8.314462681e-3 # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = 298

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 0.1
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    block = 'blockMetal'
  [../]

  [./creepSwell]
    type = ADIsoTropicHyperViscoSwellingCreep
    # absolute_tolerance = 1e-6
    # relative_tolerance = 1e-06
    hardening_exponent = 2.0
    rate_exponent = 0.15
    saturation_resistance = 2.0e-6
    hardening_modulus = 10.0e-6
    initial_resistance = 0.95e-6
    activation_energy = 37
    gas_constant = 8.314462681e-3 # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = 298

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 0.1
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    block = 'interLayer'
    enable = true

    omega = 2
    alpha1 = 0.0
    alpha2 = 1.0
    alpha3 = 0.0
    concentration = 0.0
    cref = 0.0

#    intBnd = blockCeramic_bottom
    intBnd = blockMetal_bottom

  [../]

[]

[Postprocessors]
  [./stress_yy]
    type = SideAverageValue
    variable = stress_yy
    boundary = blockMetal_top
  [../]
  [./matl_ts_min]
    type = MaterialTimeStepPostprocessor
    block = 'blockMetal'
  [../]
 [./contact]
    type = ContactDOFSetSize
    variable = frictionless_normal_lm
    subdomain = frictionless_secondary_subdomain
  [../]

[]

[VectorPostprocessors]
  [nodal_disp]
    type = NodalValueSampler
    variable = 'disp_x disp_y'
    block    = interLayer
    sort_by  = x
  []
#
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
  l_max_its = 500
  nl_max_its = 40

## Timestep setting
  dt = 0.005
  dtmin = 1e-5
  dtmax = 0.5
  end_time = 0.005
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
  execute_on = 'timestep_end'
  exodus = true
  csv    = true
  vtk = true
  file_base = rst/${name}
  [out]
    type = Checkpoint
    interval = 2
    num_files = 2
  []
  # [vtk]
  #   type = VTK
  # []
[]

[Debug]
  show_var_residual_norms = true
[]
