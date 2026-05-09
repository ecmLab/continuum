[Mesh]
  patch_size = 100
  patch_update_strategy = iteration
  [./mesh]
    type = FileMeshGenerator
    file = test_mesh_filet_fine.unv
  [../]

  [./boundary]
    type = RenameBoundaryGenerator
    input = mesh
    old_boundary_name = 'Li_bottom ceramic_top'
    new_boundary_id = '15 16'
  [../]

  # [./master]
  #   type = LowerDBlockFromSidesetGenerator
  #   input = boundary
  #   sidesets = '16'
  #   new_block_id = '3'
  # [../]
  #
  # [./slave]
  #   type = LowerDBlockFromSidesetGenerator
  #   input = boundary
  #   sidesets = '15'
  #   new_block_id = '4'
  # [../]

[]

[ICs]
  [./disp_y]
    block = 'Li_metal'
    variable = uy
    value = -2e-2
    type = ConstantIC
  [../]
[]

[GlobalParams]
  displacements = 'ux uy'
[]

[Variables]
  [./ux]
    block ='Ceramic Li_metal'
  [../]
  [./uy]
    block ='Ceramic Li_metal'
  [../]
  # [./normal_lm]
  #   block = 4
  # [../]
  # [./tangential_lm]
  #   block = 4
  #   family = MONOMIAL
  #   order = CONSTANT
  # [../]
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
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress'
    use_automatic_differentiation = true
    # use_finite_deform_jacobian = true
    block = 'Ceramic Li_metal'
  [../]
[]
[Contact]
  [./mech]
    displacements = 'ux uy'
    master = 16
    slave = 15
    system = Constraint
    formulation = mortar
    mesh = boundary
  [../]
[]
# [Constraints]
#   [./normal_lm]
#     type = NormalNodalLMMechanicalContact
#     slave = 15
#     master = 16
#     variable = normal_lm
#     master_variable = uy
#     disp_y = uy
#     ncp_function_type = min
#   [../]
#   [./normal_x]
#     type = NormalMortarMechanicalContact
#     master_boundary = 16
#     slave_boundary = 15
#     master_subdomain = 3
#     slave_subdomain = 4
#     variable = normal_lm
#     slave_variable = ux
#     component = y
#     use_displaced_mesh = true
#     compute_lm_residuals = false
#   [../]
#   [./normal_y]
#     type = NormalMortarMechanicalContact
#     master_boundary = 16
#     slave_boundary = 15
#     master_subdomain = 3
#     slave_subdomain = 4
#     variable = normal_lm
#     slave_variable = uy
#     component = y
#     use_displaced_mesh = true
#     compute_lm_residuals = false
#   [../]
#   [./tangential_x]
#     type = TangentialMortarMechanicalContact
#     master_boundary = 16
#     slave_boundary = 15
#     master_subdomain = 3
#     slave_subdomain = 4
#     variable = tangential_lm
#     slave_variable = ux
#     component = x
#     use_displaced_mesh = true
#     compute_lm_residuals = false
#   [../]
#   [./tangential_y]
#     type = TangentialMortarMechanicalContact
#     master_boundary = 16
#     slave_boundary = 15
#     master_subdomain = 3
#     slave_subdomain = 4
#     variable = tangential_lm
#     slave_variable = uy
#     component = y
#     use_displaced_mesh = true
#     compute_lm_residuals = false
#   [../]
# []

[Materials]

  [./elasticity_tensor_Li]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 7.81e-3
    poissons_ratio = 0.38
    block = 'Li_metal'
  [../]
  [./elasticity_tensor_Llzo]
    type = ADComputeIsotropicElasticityTensor
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
  automatic_scaling = true

  solve_type = 'PJFNK'
  petsc_options = '-snes_converged_reason -ksp_converged_reason -pc_svd_monitor -snes_linesearch_monitor'
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -pc_mat_solver_package'
  petsc_options_value = 'lu       NONZERO               1e-15                   superlu_dist'
  l_max_its = 30
  nl_max_its = 20


  line_search = 'none'


  nl_abs_tol = 1e-10
  # nl_rel_tol = 1e-6

  # l_max_its = 100


  start_time = 0.0
  dt = 0.001
  dtmax = 0.05
  dtmin = 1e-5


  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.001
    growth_factor = 1.5
    cutback_factor = 0.25
    optimal_iterations = 20
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
