#2D axisymmetric RZ simulation of uniaxial tension plasticity
## Using the following unit system: Length: mm, Stress: MPa

# Following the sample size used in LePage's paper
sz_x  = 1.5e-3        # axisymmetric, total width is 3mm
sz_y  = 7.5e-3        # plan symmetric, total length is 15mm
strn_rt = 2e-2     # Strain rate

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

#[Problem]
#  coord_type = RZ
#[]

[Mesh]
  type = GeneratedMesh
  xmax = ${sz_x}
  ymax = ${sz_y}
  nx = 10
  ny = 50
  dim = 2
  elem_type = QUAD4
#  second_order = true
[]

[Functions]
  [./vel2]
    type = ParsedFunction
     value  = '${sz_y} * ${strn_rt} * exp(${strn_rt}*t)'
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    volumetric_locking_correction = true
    generate_output = 'stress_yy strain_yy vonmises_stress'
    use_automatic_differentiation = true
  [../]
[]

[AuxVariables]
  [./effPlastic_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hard]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./defGrad_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]

[]

[AuxKernels]
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
  [./hard]
    type = ADMaterialRealAux
    variable = hard
    property = hardening_variable
  [../]

[]

[BCs]
  [./boty]
    type = ADDirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
    preset = true
  [../]
  [./topy]
    type = PresetVelocity
    variable = disp_y
    boundary = 'top'
    function = vel2
  [../]
  [./ux]
    type = ADDirichletBC
    variable = disp_x
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
  [./stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  [../]
  [./strain_yy]
    type = ElementAverageValue
    variable = strain_yy
  [../]
  [./u_z]
    type = AverageNodalVariableValue
    variable = disp_y
    boundary = top
  [../]
  [./effPlastic_strain]
    type = ElementAverageValue
    variable = effPlastic_strain
  [../]
  [./matl_ts_min]
    type = MaterialTimeStepPostprocessor
  [../]
  [./hard]
    type = ElementAverageValue
    variable = hard
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
  automatic_scaling = true
# Timestep setting
  dt    = 0.01
  dtmin = 1.0e-5
  dtmax = 0.5
# Solver setting
#  solve_type = 'PJFNK'
#  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart'
#  petsc_options_value = 'asm lu 1 101'
  solve_type = 'NEWTON'
  petsc_options_iname = -pc_type
  petsc_options_value = lu

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 99
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    growth_factor = 1.5
    cutback_factor = 0.85
#   optimal_iterations = 100
    timestep_limiting_postprocessor = matl_ts_min
  [../]
  end_time = 10.0
[]

[Outputs]
  exodus = true
  csv    = true
  file_base = rst/tension_constTrueStrainRate000
[]

#[Debug]
#  show_material_props = true
#[]
