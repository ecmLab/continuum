[Mesh]
  [./mesh]
    type = FileMeshGenerator
    file = test_mesh_filet_fine.unv
  [../]
  # [./corner_node]
  #   type = BoundingBoxNodeSetGenerator
  #   input = mesh
  #   new_boundary = 'monitor'
  #   top_right = '-3.1 4.9 0.0'
  #   bottom_left = '-2.9 5.1 0.0'
  # [../]
  # [./boundary]
  #   type = RenameBoundaryGenerator
  #   input = mesh
  #   old_boundary_name = 'Li_bottom'
  #   new_boundary_id = '15'
  # [../]
  # [./lower_d_block]
  #   type = LowerDBlockFromSidesetGenerator
  #   input = mesh
  #   sidesets = '15'
  #   new_block_id = '15'
  #   new_block_name = 'slave'
  # [../]
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
[Functions]
  [./pressure]
    type = PiecewiseLinear
    x = '0 0.25 1.0 5.0'
    y = '0.1 1.0 2.0 5.0'
  [../]
  [./uy]
    type = ParsedFunction
    value = '1.0*t'
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    use_displaced_mesh = true
    # volumetric_locking_correction = true
    generate_output = 'stress_yy strain_yy vonmises_stress'
    use_automatic_differentiation = true
    # use_finite_deform_jacobian = true
  [../]
[]
[Contact]
  [./mech_contact]
    disp_x = ux
    disp_y = uy
    master = 'ceramic_top'
    slave = 'Li_bottom'
    penalty = 1e-2
  [../]
[]

[Materials]

  [./elasticity_tensor_Li]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 7.81e-3
    poissons_ratio = 0.38
    block = 'Li_metal'
  [../]
  [./elasticity_tensor_Llzo]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 170e-3
    poissons_ratio = 0.3
    block = 'Ceramic'
  [../]
  [./stress_llzo]
    type = ADComputeFiniteStrainElasticStress
    block = 'Ceramic'
  [../]

  [./stress_Li]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    # perform_finite_strain_rotations = true
    block = 'Li_metal'
  [../]
  [./plas]
    type = ADIsoTropicHyperViscoCreep
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
    temperature = 398

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 1.0
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    block = 'Li_metal'
  [../]
[]

[BCs]
  [./Li_top]
    type = ADPressure
    boundary = 'Li_top'
    variable = uy
    component = 1
    function = pressure
    constant = 1.0e-6
  [../]

  [./left_right]
    type = ADPresetBC
    variable = ux
    boundary = 'ceramic_left ceramic_bot ceramic_right Li_right Li_left'
    value = 0
  [../]
  [./bottom]
    type = ADPresetBC
    variable = uy
    boundary = 'ceramic_bot'
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
    block = Li_metal
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


  line_search = 'none'


  nl_abs_tol = 5e-8
  nl_rel_tol = 1e-3

  l_max_its = 200


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
[] # Outputs
