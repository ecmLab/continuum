[Mesh]
  type = GeneratedMesh
  xmax = 1e-3
  ymax = 1e-3
  zmax = 1e-3
  nx = 1
  ny = 1
  nz = 1
  dim = 3
#  elem_type = HEX8
[]
[GlobalParams]
  displacements = 'ux uy uz'
[]

[Variables]
  [./ux]
  [../]
  [./uy]
  [../]
  [./uz]
  [../]
[]

[Functions]
  [./vel2]
    type = ParsedFunction
    # value = '1.0e-3 * (exp(3.0e-3*t) - 1.0)'
    # value  = '1.0e-3 * 2.0e-3*t'
    value = '1.0e-3 * 3e-3 * exp(3e-3*t)'
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    # volumetric_locking_correction = true
    generate_output = 'stress_zz strain_zz vonmises_stress'
    use_automatic_differentiation = true
    # use_finite_deform_jacobian = true
  [../]
[]

[BCs]
  [./symmy]
    type = ADDirichletBC
    variable = uy
    boundary = 'bottom'
    value = 0
  [../]
  [./symmx]
    type = ADDirichletBC
    variable = ux
    boundary = 'left'
    value = 0
  [../]
  [./symmz]
    type = ADDirichletBC
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
    temperature = 298

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 0.5
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
#  solve_type = 'NEWTON'
#  petsc_options_iname = -pc_type
#  petsc_options_value = lu
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -sub_pc_factor_shift_type'
  petsc_options_value = 'asm lu 2 nonzero'
  dtmax = 10.0
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  dtmin = 1.0e-4
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    growth_factor = 1.5
    cutback_factor = 0.25
    optimal_iterations = 10
    # timestep_limiting_postprocessor = matl_ts_min
  [../]
  #end_time = 333.333
   end_time = 20.0
  # num_steps = 1
[]

[Outputs]
  exodus = true
  csv = true
  file_base = rst/oneElement
[]

