[Mesh]
  [./mesh]
    type = GeneratedMeshGenerator
    xmax = 25
    ymax = 2
    nx = 50
    ny = 10
    dim = 2
    elem_type = QUAD4
  [../]
  [./m1]
    type = RenameBlockGenerator
    input = mesh
    old_block_id = '0'
    new_block_id = '1'
  [../]
  [./mesh1]
    type = GeneratedMeshGenerator
    xmax = 25
    ymax = 20
    ymin = 2
    nx = 50
    ny = 10
    dim = 2
    elem_type = QUAD4
  [../]
  [./m2]
    type = RenameBlockGenerator
    input = mesh1
    old_block_id = '0'
    new_block_id = '2'
  [../]
  [./m3]
    type = StitchedMeshGenerator
    inputs = 'm1 m2'
    stitch_boundaries_pairs = 'top bottom'
  [../]
[]


[GlobalParams]
  displacements = 'ux uy'
[]

[Variables]
  [./ux]
  [../]
  [./uy]
  [../]
[]


[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    use_displaced_mesh = true
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy vonmises_stress strain_xx strain_yy elastic_strain_xx elastic_strain_yy'
    # use_finite_deform_jacobian = true
    use_automatic_differentiation = true
  [../]
[]



[AuxVariables]
  [./conc]
  [../]

  [./Fp_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plastic_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strength]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./conc]
    type = FunctionAux
    variable = conc
    function = '1.0e-15*t'
  [../]

  [./fp_yy]
    type = ADRankTwoAux
    variable = Fp_yy
    rank_two_tensor = plastic_distortion
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./peeq]
    type = ADMaterialRealAux
    variable = plastic_strain
    property = effective_plastic_strain
    execute_on = timestep_end
  [../]
  [./strength]
    type = ADMaterialRealAux
    variable = strength
    property = yield_strength
    execute_on = 'TIMESTEP_END'
  [../]
[]

[BCs]
  [./symmy]
    type = ADDirichletBC
    variable = uy
    boundary = bottom
    value = 0
    preset = true
  [../]
  # [./topy]
  #   type = PresetBC
  #   variable = uy
  #   boundary = top
  #   value = 0
  # [../]
  [./symmx]
    type = ADDirichletBC
    variable = ux
    boundary = left
    value = 0
    preset = true
  [../]


  # [./tdisp]
  #   type = FunctionPresetBC
  #   variable = uy
  #   boundary = top
  #   function = uy
  # [../]
[]


[Materials]

  [./elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 5e-3
    poissons_ratio = 0.3
  [../]
  [./stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    block = '1'
  [../]
  [./plas]
    type = ADIsoTropicHyperViscoSwelling
    absolute_tolerance = 1e-5
    hardening_exponent = 1.8
    saturation_resistance = 8.0e-6
    initial_resistance = 2.0e-6
    hardening_modulus = 40.0e-6
    rate_exponent = 0.18
    reference_strain_rate = 0.05
    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 1.0
    omega = 1e13
    alpha = '0 1 0'
    concentration = conc
    cref = 0.0
    block = 1
    # internal_solve_full_iteration_history = true
    # internal_solve_output_on = always
  [../]
  [./stress2]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas2'
    block = '2'
  [../]
  [./plas2]
    type = ADIsoTropicHyperViscoSwelling
    absolute_tolerance = 1e-5
    hardening_exponent = 1.8
    saturation_resistance = 8.0e-6
    initial_resistance = 2.0e-6
    hardening_modulus = 40.0e-6
    rate_exponent = 0.18
    reference_strain_rate = 0.05
    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 1.0
    omega = 1e-3
    alpha = '0 1 0'
    concentration = conc
    cref = 0.0
    block = 2
    # internal_solve_full_iteration_history = true
    # internal_solve_output_on = always
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
  [./u_y]
    type = AverageNodalVariableValue
    variable = uy
    boundary = top
  [../]
  [./peeq]
    type = ElementAverageValue
    variable = plastic_strain
  [../]
  [./fp_yy]
    type = ElementAverageValue
    variable = Fp_yy
  [../]
  [./stregnth]
    type = ElementAverageValue
    variable = strength
  [../]
  [./elastic_strain_yy]
    type = ElementAverageValue
    variable = elastic_strain_yy
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient

  #Preconditioned JFNK (default)
  automatic_scaling = true
  dt = 1.0
  solve_type = 'NEWTON'
  petsc_options_iname = -pc_type
  petsc_options_value = lu
  dtmax = 5
  nl_rel_tol = 1e-3
  nl_abs_tol = 1e-10
  dtmin = 1.0e-8
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    growth_factor = 1.5
    cutback_factor = 0.25
    optimal_iterations = 25
  [../]
  # num_steps = 10
  end_time = 200.0
[]

[Outputs]
  exodus = true
  csv = true
[]

[Debug]
  show_material_props = true
[]
