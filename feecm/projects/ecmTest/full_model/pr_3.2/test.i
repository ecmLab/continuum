strnRt = 1e-5

[Mesh]
  patch_size = 80
  patch_update_strategy = auto
  parallel_type = REPLICATED
  [./blocks]
    type = FileMeshGenerator
    file = data/cosine_conformal_full.msh
  [../]

[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [./disp_x]
  [../]

  [./disp_y]
  [../]

[]

[AuxVariables]
#  [./Fp_yy]
#    order = CONSTANT
#    family = MONOMIAL
#    block = 'blockMetal interLayer'
#  [../]
  [./plastic_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_rate]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strength]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]
  [./pressure]
    # type = ParsedFunction
    # value = '0.1e-6'
    type = PiecewiseLinear
    data_file = 'pressure_time.csv'
    format = columns
  [../]

  [./uy]
    type = ParsedFunction
    value = '0.01*t'
  [../]
  [./vel]
    type = ParsedFunction
    value = '${strnRt}*exp(${strnRt}*t)'
  [../]

[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'disp_x disp_y mech_contact_normal_lm'
  acceptable_iterations = 2
[]

[Modules/TensorMechanics/Master]
  [./action]
    add_variables = true
    strain = FINITE
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy vonmises_stress strain_xx strain_yy'
    block = 'blockCeramic blockMetal interLayer'
    extra_vector_tags = 'ref'
    use_automatic_differentiation = true
  [../]
[]

# [Problem]
#   type = ReferenceResidualProblem
#   # solution_variables = 'ux uy normal_lm thermal_lm li_ion_V thermal_lm2 li_metal_conc'
#   extra_tag_vectors = 'ref'
#   reference_vector = 'ref'
#   group_variables = 'ux uy mech_contact_normal_lm'
#   acceptable_iterations = 2
# []

[Contact]
  [./mech_contact]
    mesh = blocks
    secondary = blockMetal_bottom
    primary = blockCeramic_top
    # penalty = 1e-2
#    formulation = kinematic
#    normal_smoothing_distance = 0.5
#    tangential_tolerance = 0.25
    formulation = mortar
  [../]
[]

[AuxKernels]

#  [./fp_yy]
#    type = ADRankTwoAux
#    variable = Fp_yy
#    rank_two_tensor = plastic_distortion
#    index_i = 1
#    index_j = 1
#    execute_on = timestep_end
#     block = 'interLayer blockMetal'
#  [../]
  [./peeq]
    type = ADMaterialRealAux
    variable = plastic_strain
    property = effective_plastic_strain
    execute_on = timestep_end
     block = 'interLayer blockMetal'
  [../]
  [./strain_rate]
    type = ADMaterialRealAux
    variable = strain_rate
    property = plastic_strain_rate
    block = 'interLayer blockMetal'
  [../]
  [./strength]
    type = ADMaterialRealAux
    variable = strength
    property = yield_strength
    execute_on = 'TIMESTEP_END'
    block = 'interLayer blockMetal'
  [../]
[]

[BCs]
## Mechanics
#  [./pressure]
#    type = ADPressure
#    variable = disp_y
#    component = 1
#    function = pressure
#    boundary = 'blockCeramic_bottom'
#    use_displaced_mesh = true
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
    function = vel
  [../]
  [./left_right]
    type = ADDirichletBC
    variable = disp_x
    boundary = 'blockCeramic_left blockCeramic_bottom blockCeramic_right blockMetal_right blockMetal_left'
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
    block = 'blockMetal interLayer'
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
    block = 'blockMetal interLayer'
  [../]

[]

[Postprocessors]
  [./stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  [../]
  [./strain_yy]
    type = ElementAverageValue
    variable = strain_yy
  [../]
  [./peeq]
    type = ElementAverageValue
    variable = plastic_strain
  [../]
#  [./fp_yy]
#    type = ElementAverageValue
#    variable = Fp_yy
#  [../]
  [./stregnth]
    type = ElementAverageValue
    variable = strength
  [../]
  [./ext_pressure]
    type = SideAverageValue
    variable = stress_yy
    boundary = 'blockMetal_top'
  [../]
  [./matl_ts_min]
    type = MaterialTimeStepPostprocessor
    block = 'interLayer blockMetal'
  [../]
[]

[VectorPostprocessors]
  [nodal_disp]
    type = NodalValueSampler
    variable = 'disp_x disp_y'
    block    = interLayer
    sort_by  = x
  []

[]

[Executioner]
  type = Transient
## Solver setting
  solve_type = 'NEWTON'
  automatic_scaling = true
#  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor -snes_ksp_ew'
#  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_force_iteration -pc_factor_shift_type -pc_factor_shift_amount -mat_mffd_err'
#  petsc_options_value = 'lu superlu_dist 1     NONZERO               1e-15               1e-5'
  petsc_options = '-snes_converged_reason -ksp_converged_reason'
  petsc_options_iname = '-pc_type -mat_mffd_err -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu       1e-5          NONZERO               1e-15'

#  solve_type = 'PJFNK'
#  automatic_scaling = true
#  petsc_options = '-snes_ksp_ew'
#  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
#  petsc_options_value = 'hypre    boomeramg      101'

  line_search = 'none'
  nl_abs_tol = 1e-15
  nl_rel_tol = 1e-9
  l_max_its = 100
  nl_max_its = 35

## Timestep setting
  dt = 0.0001
  dtmin = 1e-5
  dtmax = 0.25
  end_time = 100.0
  timestep_tolerance = 1e-6
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.0001
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 30
    timestep_limiting_postprocessor = matl_ts_min
  [../]

[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
  csv    = true
  file_base = cosine_conformal_full
  [out]
    type = Checkpoint
    interval = 2
    num_files = 2
  []
[]

[Debug]
  show_var_residual_norms = true
  show_material_props = true
[]
