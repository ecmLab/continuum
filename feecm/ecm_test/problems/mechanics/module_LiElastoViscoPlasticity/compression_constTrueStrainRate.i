#2D axisymmetric RZ simulation of uniaxial tension plasticity

[GlobalParams]
  displacements = 'disp_r disp_z'
[]

[Problem]
  coord_type = RZ
[]

[Mesh]
  type = GeneratedMesh
  xmax = 6.35e-3
  ymax = 15.5e-3
  nx = 15
  ny = 25
  dim = 2
  second_order = true
[]

[Functions]
  [./vel2]
    type = ParsedFunction
     value  = '-1.0e-2 * 3.0e-3 * exp(3.0e-3*t)'
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    volumetric_locking_correction = true
    generate_output = 'stress_yy strain_yy vonmises_stress'
    use_automatic_differentiation = true
    use_finite_deform_jacobian = true
  [../]
[]

[AuxVariables]
  [./Fp_yy]
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
  [./mandel_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./logStrain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./defGrad_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]

[]

[AuxKernels]
  [./fp_yy]
    type = ADRankTwoAux
    variable = Fp_yy
    rank_two_tensor = plastic_distortion
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./log_strain_yy]
    type = ADRankTwoAux
    variable = logStrain_yy
    rank_two_tensor = logarithmic_elastic_strain
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./def_grad_yy]
    type = ADRankTwoAux
    variable = defGrad_yy
    rank_two_tensor = deformation_gradient
    index_i = 1
    index_j = 1
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
#  [./effective_strain]
#    type = ADRankTwoScalarAux
#    rank_two_tensor = total_strain
#    variable = effstrain
#    scalar_type = EffectiveStrain
#  [../]
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
  [./mandel_yy]
    type = ADRankTwoAux
    variable = mandel_yy
    rank_two_tensor = mandel_stress
    index_i = 1
    index_j = 1
  [../]

[]

[BCs]
  [./boty]
    type = ADDirichletBC
    variable = disp_z
    boundary = 'bottom'
    value = 0
    preset = true
  [../]
  [./topy]
    type = PresetVelocity
    variable = disp_z
    boundary = 'top'
    function = vel2
  [../]
  [./ux]
    type = ADDirichletBC
    variable = disp_r
    boundary = left
    value = 0
    preset = true
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 7.81e3
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
    saturation_resistance = 2.0
    initial_resistance = 0.95
    hardening_modulus = 10.0
    rate_exponent = 0.15
    activation_energy = 37
    gas_constant = 8.314462681e-3 # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = 298

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 0.5
#    internal_solve_full_iteration_history = true
#    internal_solve_output_on = on_error
  [../]
[]

[Postprocessors]
#  [./stress_yy]
#    type = SideAverageValue
#    variable = stress_yy
#    boundary = top
#  [../]
#  [./strain_zz]
#    type = SideAverageValue
#    variable = strain_zz
#    boundary = top
#  [../]
  [./effPlastic_strain]
    type = ElementAverageValue
    variable = effPlastic_strain
  [../]
  [./yield_stregnth]
    type = ElementAverageValue
    variable = yldStrength
  [../]
  [./matl_ts_min]
    type = MaterialTimeStepPostprocessor
  [../]
#  [./eff_strain]
#    type=ElementAverageValue
#    variable = effstrain
#  [../]
  [./von_mises]
    type = SideAverageValue
    variable = vonmises_stress
    boundary = top
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
  automatic_scaling = true
# Timestep setting
  dt    = 0.01
  dtmin = 1.0e-8
  dtmax = 0.1
  end_time = 10.0
# Solver setting
#  solve_type = 'PJFNK'
#  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart'
#  petsc_options_value = 'asm lu 1 101'
  solve_type = 'NEWTON'
  petsc_options_iname = -pc_type
  petsc_options_value = lu

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  nl_max_its = 15
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    growth_factor = 1.5
    cutback_factor = 0.25
    optimal_iterations = 15
    timestep_limiting_postprocessor = matl_ts_min
  [../]
[]

[Outputs]
  exodus = true
  csv    = true
  file_base = rst/compression_constTrueStrainRate
[]

#[Debug]
#  show_material_props = true
#[]
