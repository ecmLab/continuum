[Mesh]
  type = GeneratedMesh
  xmax = 1e-3
  ymax = 1e-3
  zmax = 1e-3
  nx = 1
  ny = 1
  nz = 1
  dim = 3
  elem_type = HEX8
[]
[GlobalParams]
  displacements = 'ux uy uz'
[]

[Variables]
  [./ux]
    # [./InitialCondition]
    #   type = RandomIC
    #   max = 1.0e-6
    #   min = -1.0e-6
    # [../]
  [../]
  [./uy]
    # [./InitialCondition]
    #   type = RandomIC
    #   max = 1.0e-6
    #   min = -1.0e-6
    # [../]
  [../]
  [./uz]
    # [./InitialCondition]
    #   type = RandomIC
    #   max = 1.0e-6
    #   min = -1.0e-6
    # [../]
  [../]
[]

[Functions]
  [./vel2]
    type = ParsedFunction
    # value = '1.0e-3 * (exp(3.0e-3*t) - 1.0)'
    # value  = '1.0e-3 * 2.0e-3*t'
    value = '1.0e-3 * 2e-2 * exp(2e-2*t)'
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    use_displaced_mesh = true
    # volumetric_locking_correction = true
    generate_output = 'stress_zz strain_zz vonmises_stress'
    use_automatic_differentiation = true
    # use_finite_deform_jacobian = true
  [../]
[]

[AuxVariables]
  [./Fp_zz]
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
  [./effstrain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_rate]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hard]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./mandel]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./log_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./def_grad]
    order = CONSTANT
    family = MONOMIAL
  [../]


[]

[AuxKernels]
  [./fp_zz]
    type = RankTwoAux
    variable = Fp_zz
    rank_two_tensor = plastic_distortion
    index_i = 2
    index_j = 2
    execute_on = timestep_end
  [../]
  [./log_strain_zz]
    type = RankTwoAux
    variable = log_strain
    rank_two_tensor = logarithmic_elastic_strain
    index_i = 2
    index_j = 2
    execute_on = timestep_end
  [../]

  [./def_grad_zz]
    type = RankTwoAux
    variable = def_grad
    rank_two_tensor = deformation_gradient
    index_i = 2
    index_j = 2
    execute_on = timestep_end
  [../]

  [./peeq]
    type = MaterialRealAux
    variable = plastic_strain
    property = effective_plastic_strain
    execute_on = timestep_end
  [../]
  [./strength]
    type = MaterialRealAux
    variable = strength
    property = yield_strength
    execute_on = 'TIMESTEP_END'
  [../]
  [./effective_strain]
    type = RankTwoScalarAux
    rank_two_tensor = total_strain
    variable = effstrain
    scalar_type = EffectiveStrain
  [../]

  [./strain_rate]
    type = MaterialRealAux
    variable = strain_rate
    property = plastic_strain_rate
  [../]
  [./hard]
    type = MaterialRealAux
    variable = hard
    property = hardening_variable
  [../]
  [./mandel_zz]
    type = RankTwoAux
    variable = mandel
    rank_two_tensor = mandel_stress
    index_i = 2
    index_j = 2
  [../]


[]

[BCs]
  [./symmy]
    type = ADPresetBC
    variable = uy
    boundary = 'bottom'
    value = 0
  [../]
  [./symmx]
    type = ADPresetBC
    variable = ux
    boundary = 'left'
    value = 0
  [../]
  [./symmz]
    type = ADPresetBC
    variable = uz
    boundary = back
    value = 0
  [../]
  [./tdisp]
    type = PresetVelocity
    variable = uz
    boundary = front
    function = vel2
  [../]
[]


[Materials]

  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 7.81e9
    poissons_ratio = 0.38
  [../]
  [./stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    perform_finite_strain_rotations = true
  [../]
  [./plas]
    type = ADIsoTropicHyperViscoCreep
    absolute_tolerance = 1e-8
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
    temperature = 198

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 1.0
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
  [../]
[]

[Postprocessors]
  [./stress_zz]
    type = ElementAverageValue
    variable = stress_zz
  [../]
  [./strain_zz]
    type = ElementAverageValue
    variable = strain_zz
  [../]
  [./u_z]
    type = AverageNodalVariableValue
    variable = uz
    boundary = front
  [../]
  [./peeq]
    type = ElementAverageValue
    variable = plastic_strain
  [../]
  [./fp_zz]
    type = ElementAverageValue
    variable = Fp_zz
  [../]
  [./stregnth]
    type = ElementAverageValue
    variable = strength
  [../]
  [./matl_ts_min]
    type = MaterialTimeStepPostprocessor
  [../]
  [./eff_strain]
    type=ElementAverageValue
    variable = effstrain
  [../]

  [./pl_rate]
    type = ElementAverageValue
    variable = strain_rate
  [../]

  [./hard]
    type = ElementAverageValue
    variable = hard
  [../]
  [./mandel]
    type = ElementAverageValue
    variable = mandel
  [../]
  [./log_elastic]
    type = ElementAverageValue
    variable = log_strain
  [../]

  [./def_grad]
    type = ElementAverageValue
    variable = def_grad
  [../]
  [./von_mises]
    type = ElementAverageValue
    variable = vonmises_stress
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
  dt = 0.01
  solve_type = 'NEWTON'
  petsc_options_iname = -pc_type
  petsc_options_value = lu
  dtmax = 0.1
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  dtmin = 1.0e-8
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    growth_factor = 1.5
    cutback_factor = 0.25
    optimal_iterations = 10
    # timestep_limiting_postprocessor = matl_ts_min
  [../]
  end_time = 50.0
  # end_time = 20.0
  # num_steps = 1
[]

[Outputs]
  exodus = true
  csv = true
[]

[Debug]
  show_material_props = true
[]
