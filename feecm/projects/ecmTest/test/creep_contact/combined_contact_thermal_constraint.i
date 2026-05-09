[Mesh]
  dim = 2
  displacements = 'ux uy'
  patch_size = 500
  patch_update_strategy = iteration
  [./mesh]
    type = FileMeshGenerator
    file = test_mesh_filet_fine.unv
  [../]
  [./slave_boundary]
    type = RenameBoundaryGenerator
    input = mesh
    old_boundary_name = 'Li_bottom'
    new_boundary_id = '15'
  [../]
  [./master_boundary]
    type = RenameBoundaryGenerator
    input = slave_boundary
    old_boundary_name = 'ceramic_top'
    new_boundary_id = '16'
  [../]
  [./slave_block]
    type = LowerDBlockFromSidesetGenerator
    input = master_boundary
    sidesets = '15'
    new_block_id = '3'
    new_block_name = 'slave'
  [../]
  [./master_block]
    type = LowerDBlockFromSidesetGenerator
    input = slave_block
    sidesets = '16'
    new_block_id = '4'
    new_block_name = 'master'
  [../]
[]

[GlobalParams]
  displacements = 'ux uy'
[]

[Variables]
  [./ux]
    block = 'Ceramic Li_metal'
  [../]
  [./uy]
    block = 'Ceramic Li_metal'
  [../]
  [./temp]
    block = 'Ceramic Li_metal'
  [../]
  [./thermal_lm]
    block = '4'
  [../]
[]
[Functions]
  [./pressure]
    type = ConstantFunction
    value = 1.0e-6
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
    generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress'
    use_automatic_differentiation = true
    block = 'Ceramic Li_metal'
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
[Constraints]
  [./thermal_constraint]
    type = GapDisplacementConductanceConstraint
    variable = thermal_lm
    slave_variable = temp
    slave_boundary = 16
    slave_subdomain = 4
    master_boundary = 15
    master_subdomain = 3
    gap_limit = 1.0e-6
    k = 100.0
  [../]
[]
[Kernels]
  [./temp]
    type = ADHeatConduction
    diffusion_coefficient = thermal_conductivity
    block = 'Ceramic Li_metal'
    variable = temp
  [../]
  [./temp_dt]
    type = ADHeatConductionTimeDerivative
    density_name = density
    variable = temp
    block = 'Ceramic Li_metal'
    specific_heat = 1.0
  [../]
[]
[Materials]
  [./thermal_conductivity]
    type = HeatConductionMaterial
    thermal_conductivity = 1.0
    block = 'Ceramic Li_metal'
  [../]
  [./density]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = '1.0'
    block = 'Ceramic Li_metal'
  [../]

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
    temperature = 298

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 1.0
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    block = 'Li_metal'
  [../]
[]

[BCs]
  [./Li_top]
    type = ADFunctionPresetBC
    boundary = 'ceramic_bot'
    variable = uy
    component = 1
    function = uy
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
    boundary = 'Li_top'
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
  [./bot_stress]
    type = SideAverageValue
    variable = stress_yy
    boundary = 'ceramic_bot'
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


  nl_abs_tol = 5e-8
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
