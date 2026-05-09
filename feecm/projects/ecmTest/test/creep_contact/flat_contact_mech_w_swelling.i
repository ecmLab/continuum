
[Mesh]
  [./mesh1]
    type = GeneratedMeshGenerator
    xmax = 1e-3
    xmin = -1.0e-3
    ymax = 0.5e-3
    ymin = 0.0
    dim = 2
    nx = 25
    ny = 15
  [../]
  [./m1]
    type = RenameBlockGenerator
    input = mesh1
    old_block_id = '0'
    new_block_id = '1'
  [../]
  [./slave_boundary]
    type = RenameBoundaryGenerator
    input = m1
    old_boundary_name = 'bottom right top left'
    new_boundary_id = '16 17 18 19'
  [../]

  [./mesh2]
    type = GeneratedMeshGenerator
    xmax = 1e-3
    xmin = -1.0e-3
    ymax = 0.0
    ymin = -0.5e-3
    dim = 2
    nx = 15
    ny = 15
  [../]

  [./m2]
    type = RenameBlockGenerator
    input = mesh2
    old_block_id = '0'
    new_block_id = '2'
  [../]

  [./master_boundary]
    type = RenameBoundaryGenerator
    input = m2
    old_boundary_name = 'bottom right top left'
    new_boundary_id = '11 12 15 14'
  [../]
  [./cmg]
    type = CombinerGenerator
    inputs = 'master_boundary slave_boundary'
  [../]

[]


[GlobalParams]
  displacements = 'ux uy'
[]

[Variables]
  [./ux]
    block = '1 2'
  [../]
  [./uy]
    block = '1 2'
  [../]
[]

[AuxVariables]
  [./conc]
    block = '2'
  [../]
[]
[AuxKernels]
  [./conc]
    type = FunctionAux
    variable = conc
    function = '1.0e-15*t'
    block = '2'
  [../]
[]

[Functions]
  [./pressure]
    type = ConstantFunction
    value = 1.0e-6
  [../]
  [./uy]
    type = ConstantFunction
    value = '0.0'
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    use_displaced_mesh = true
    # volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress'
    use_automatic_differentiation = true
    block = '1 2'
    # use_finite_deform_jacobian = true
  [../]
[]
[Contact]
  [./mech_contact]
    disp_x = ux
    disp_y = uy
    master = 15
    slave = 16
    penalty = 1e-2
  [../]
[]
[Materials]
  [./elasticity_tensor_Li]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 7.81e-3
    poissons_ratio = 0.38
    block = '2'
  [../]
  [./elasticity_tensor_Llzo]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 170e-3
    poissons_ratio = 0.3
    block = '1'
  [../]
  [./stress_llzo]
    type = ADComputeFiniteStrainElasticStress
    block = '1'
  [../]

  [./stress_Li]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    # perform_finite_strain_rotations = true
    block = '2'
  [../]
  [./plas]
    type = ADIsoTropicHyperViscoSwellingCreep
    absolute_tolerance = 1e-6
    # relative_tolerance = 1e-06
    hardening_exponent = 2.0
    saturation_resistance = 2.0e-6
    initial_resistance = 0.95e-6
    hardening_modulus = 10.0e-6
    rate_exponent = 0.15
    activation_energy = 37
    gas_constant = 8.314462681e-3 # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = 298
    omega = 1e13
    alpha1 = 0.0
    alpha2 = 1.0
    alpha3 = 0.0
    concentration = conc
    cref = 0.0

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 1.0
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    block = '2'
  [../]
[]

[BCs]
  [./Li_top]
    type = ADPresetBC
    boundary = '18'
    variable = uy
    value =  0.0
    # component = 1
    # function = uy
  [../]

  [./left_right]
    type = ADPresetBC
    variable = ux
    boundary = '11 12 14 17 19'
    value = 0
  [../]
  [./bottom]
    type = ADPresetBC
    variable = uy
    boundary = '11'
    value = 0
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./matl_ts_min]
    type = MaterialTimeStepPostprocessor
    block = 2
  [../]
  [./bot_stress]
    type = SideAverageValue
    variable = stress_yy
    boundary = 18
  [../]
[]

[Executioner]

  type = Transient

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'
  automatic_scaling = true
  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre    boomeramg      101'
  # petsc_options_iname = -pc_type
  # petsc_options_value = lu


  line_search = 'none'


  nl_abs_tol = 1e-9
  # nl_rel_tol = 1e-6

  l_max_its = 100


  start_time = 0.0
  dt = 0.001
  dtmax = 0.05
  dtmin = 1e-5
  end_time = 2.0
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.001
    growth_factor = 1.5
    cutback_factor = 0.25
    optimal_iterations = 40
    timestep_limiting_postprocessor = matl_ts_min
  [../]

[] # Executioner

[Outputs]
  [./out]
    type = Exodus
    elemental_as_nodal = true
  [../]
  [./csv]
    type = CSV
  [../]
[] # Outputs
