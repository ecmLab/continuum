elm_sz  = 1e-3
strn_rt = 4e-5

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  type = GeneratedMesh
  xmax = ${elm_sz}
  ymax = ${elm_sz}
  zmax = ${elm_sz}
  nx = 1
  ny = 1
  nz = 1
  dim = 3
  elem_type = HEX8
[]

[Functions]
  [./vel2]
    type = ParsedFunction
    value = '${elm_sz} * ${strn_rt} * exp(${strn_rt}*t)'
  [../]
[]

[Variables]
  [./disp_x]
  [../]

  [./disp_y]
  [../]

  [./disp_z]
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
  [./effPlastic_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./yldStrength]
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
  [./mandel_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./logStrain_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./defGrad_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]

[]

[AuxKernels]
  [./fp_zz]
    type = ADRankTwoAux
    variable = Fp_zz
    rank_two_tensor = plastic_distortion
    index_i = 2
    index_j = 2
    execute_on = timestep_end
  [../]
  [./log_strain_zz]
    type = ADRankTwoAux
    variable = logStrain_zz
    rank_two_tensor = logarithmic_elastic_strain
    index_i = 2
    index_j = 2
    execute_on = timestep_end
  [../]

  [./def_grad_zz]
    type = ADRankTwoAux
    variable = defGrad_zz
    rank_two_tensor = deformation_gradient
    index_i = 2
    index_j = 2
    execute_on = timestep_end
  [../]

  [./effective_plastic_strain]
    type = ADMaterialRealAux
    variable = effPlastic_strain
    property = effective_plastic_strain
    execute_on = timestep_end
  [../]
  [./yield_strength]
    type = ADMaterialRealAux
    variable = yldStrength
    property = yield_strength
    execute_on = 'TIMESTEP_END'
  [../]
  [./effective_strain]
    type = ADRankTwoScalarAux
    rank_two_tensor = total_strain
    variable = effstrain
    scalar_type = EffectiveStrain
  [../]

  [./strain_rate]
    type = ADMaterialRealAux
    variable = strain_rate
    property = plastic_strain_rate
  [../]
  [./hard]
    type = ADMaterialRealAux
    variable = hard
    property = hardening_variable
  [../]
  [./mandel_zz]
    type = ADRankTwoAux
    variable = mandel_zz
    rank_two_tensor = mandel_stress
    index_i = 2
    index_j = 2
  [../]

[]

[BCs]
  [./symmy]
    type = ADDirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
    preset = true
  [../]
  [./symmx]
    type = ADDirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
    preset = true
  [../]
  [./symmz]
    type = ADDirichletBC
    variable = disp_z
    boundary = back
    value = 0
    preset = true
  [../]
  [./tdisp]
    type = PresetVelocity
    variable = disp_z
    boundary = front
    function = vel2
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 7.81e9
    poissons_ratio = 0.38
  [../]
  [./stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    perform_finite_strain_rotations = true
  [../]
  [./plas]
    type = ADIsoTropicHyperVisco
    absolute_tolerance = 1e-8
    # relative_tolerance = 1e-06
    hardening_exponent = 2.0
    saturation_resistance = 2.0e6
    initial_resistance = 0.95e6
    hardening_modulus = 10.0e6
    rate_exponent = 0.15
    reference_strain_rate = 0.05

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
    variable = disp_z
    boundary = front
  [../]
  [./effPlastic_strain]
    type = ElementAverageValue
    variable = effPlastic_strain
  [../]
  [./fp_zz]
    type = ElementAverageValue
    variable = Fp_zz
  [../]
  [./yield_stregnth]
    type = ElementAverageValue
    variable = yldStrength
  [../]
  [./matl_ts_min]
    type = MaterialTimeStepPostprocessor
  [../]
  [./eff_strain]
    type=ElementAverageValue
    variable = effstrain
  [../]
#  [./sat]
#    type = ElementAverageValue
#    variable = saturation
#  [../]

  [./strain_rate]
    type = ElementAverageValue
    variable = strain_rate
  [../]

  [./hard]
    type = ElementAverageValue
    variable = hard
  [../]

  [./mandel_zz]
    type = ElementAverageValue
    variable = mandel_zz
  [../]
  [./logStrain_zz]
    type = ElementAverageValue
    variable = logStrain_zz
  [../]

  [./defGrad_zz]
    type = ElementAverageValue
    variable = defGrad_zz
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
#  dtmax = 1.0
  dtmax = 75.00000000000000000000
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-10
  dtmin = 1.0e-8
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    growth_factor = 1.5
    cutback_factor = 0.25
    optimal_iterations = 10
    timestep_limiting_postprocessor = matl_ts_min
  [../]
#  end_time = 12.5
  end_time = 12500.00000000000000000000
  # num_steps = 100
[]

[Outputs]
  exodus = true
  csv = true
  file_base = rst/tension_constTrueStrainRateR4
[]

[Debug]
  show_material_props = true
[]
