# This is a model with 2 flat blocks
# block 1 - Ceramic Electrolyte (top)
# block 2 - Li_metal with Anand model swelling
# Li_ion_conc -- modeled as a steady state thermal problem
# All units in um
elem = QUAD4
order = FIRST
[Mesh]
  patch_size = 80
  patch_update_strategy = auto
  displacements = 'ux uy'
  parallel_type = REPLICATED
  [./film]
    type = GeneratedMeshGenerator
    xmax = 20
    xmin = 0
    ymax = 0.1
    ymin = 0.0
    dim = 2
    nx = 50
    ny = 4
    elem_type = ${elem}
  [../]
  [./film1]
    type = RenameBlockGenerator
    input = film
    old_block_id = '0'
    new_block_id = '5'
  [../]

  [./electrolyte1]
    type = GeneratedMeshGenerator
    xmax = 20
    xmin = 0
    ymax = 25
    ymin = 0.1
    dim = 2
    nx = 50
    ny = 15
    elem_type = ${elem}
    # bias_y = 1.5
  [../]
  [./m1]
    type = RenameBlockGenerator
    input = electrolyte1
    old_block_id = '0'
    new_block_id = '1'
  [../]

  [./electrolyte]
    type = StitchedMeshGenerator
    inputs = 'film1 m1'
    stitch_boundaries_pairs = 'top bottom'
    clear_stitched_boundary_ids = true
  [../]

  [./slave_boundary]
    type = RenameBoundaryGenerator
    input = electrolyte
    old_boundary_name = 'bottom right top left'
    new_boundary_id = '16 17 18 19'
  [../]

  [./Li_1]
    type = GeneratedMeshGenerator
    xmax = 20
    xmin = 0
    ymax = 0.0
    ymin = -0.5
    dim = 2
    nx = 60
    ny = 5
    elem_type = ${elem}
  [../]

  [./Li]
    type = RenameBlockGenerator
    input = Li_1
    old_block_id = '0'
    new_block_id = '2'
  [../]
  [./Li_2]
    type = GeneratedMeshGenerator
    xmax = 20
    xmin = 0
    ymax = -0.5
    ymin = -30
    dim = 2
    nx = 60
    ny = 15
    elem_type = ${elem}
    bias_y = 0.8333
  [../]

  [./m3]
    type = RenameBlockGenerator
    input = Li_2
    old_block_id = '0'
    new_block_id = '8'
  [../]
  [./m2]
    type = StitchedMeshGenerator
    inputs = 'Li m3'
    stitch_boundaries_pairs = 'bottom top'
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
  # [./slave_boundary_block]
  #   type = LowerDBlockFromSidesetGenerator
  #   input = cmg
  #   sidesets = '16'
  #   new_block_id = '4'
  # [../]
  # [./master_boundary_block]
  #   type = LowerDBlockFromSidesetGenerator
  #   input = slave_boundary_block
  #   sidesets = '15'
  #   new_block_id = '3'
  # [../]

[]


[GlobalParams]
  displacements = 'ux uy'
[]

[Variables]
  [./ux]
    block = '1 2 5 8'
    order = ${order}
    # scaling = 1e12
  [../]
  [./uy]
    block = '1 2 5 8'
    order = ${order}
    # scaling = 1e12
  [../]


  # [./normal_lm]
  #   block = 3
  #   order = ${order}
  # [../]


[]

[Functions]
  [./pressure]
    type = ParsedFunction
    # value = '1.2e-6'
    # type = ParsedFunction
    value = 'if (t <= 10.0, 0.04e-6*t + 0.4e-6, 1.4e-6)'
  [../]
  [./flux]
    type = ParsedFunction
    value = 'if (t > 10, if (t>= 20, -4e-12, -4e-12*(t-10)/10.0), 5e-10)'
    # value = '-2e-12/(1.0 + exp(-2.0*0.6*(t-5.0)))'
  [../]
  [./gapk]
    # type = ParsedFunction
    # value = ' 1e-9*(1.0/(1.0 + exp(2.0*3e4*(x-1e-4))))'
    type = PiecewiseLinear
    x = '-10.0 0 1e-5 1e-4 1 10 100'
    y = '1e-9 1e-9 1e-9 0 0 0 0'
  [../]
  # [./gapk]
  #   type = ParsedFunction
  #   vars = 'k1'
  #   vals = 'gapk_1'
  #   value = 'if (t <= 10, 1e-9, k1)'
  # [../]
  [./gapk1]
    # type = ParsedFunction
    # value = ' 1e3*(1.0/(1.0 + exp(2.0*3e4*(x-1e-4))))'
    type = PiecewiseLinear
    x = '-10.0 0 1e-5 1e-4 1 10 100'
    y = '1e-2 1e-2 1e-2 0 0 0 0'
  [../]
  # [./gapk1]
  #   type = ParsedFunction
  #   vars = 'k11'
  #   vals = 'gapk_11'
  #   value = 'if (t<=10, 1e1, k11)'
  # [../]

  [./k_function]
    type = ParsedFunction
    # value = 'if (abs(x) > 20.0, 1e-9, if (cos(2*3.1415*x/20.0)/2.0 < 0, 1e-11, 1e-9))'
    # value = 'if (abs(x) <= 1.5, 1e-9, 1e-11)'
    # value = '1.0e-9*exp(-x^2/0.1) + 1e-11'
    # value = '1e-9'
    # value = 'if (t <= 10, 1e-9, if (abs(x) <= 1.5, 1e-9, 1e-11))'
    value = 'if (t <=10, 1e-8,(1.0/(1.0 + exp(2.0*2.5*(x- 5.0))))*(1e-10-1e-11) + 1e-11)'

  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    use_displaced_mesh = true
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress hydrostatic_stress'
    use_automatic_differentiation = true
    block = '1 2 5 8'
    displacements = 'ux uy'
    extra_vector_tags = 'ref'
    # use_finite_deform_jacobian = true
  [../]
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'ux uy'
  acceptable_iterations = 2
[]

# [Constraints]
#   [./disp_x]
#     type = RANFSNormalMechanicalContact
#     slave = 15
#     master = 16
#     variable = ux
#     master_variable = ux
#     component = x
#   [../]
#   [./disp_y]
#     type = RANFSNormalMechanicalContact
#     slave = 15
#     master = 16
#     variable = uy
#     master_variable = uy
#     component = y
#   [../]
# []


[Contact]
  [./mech]
    master = 16
    slave = 15
    displacements = 'ux uy'
    capture_tolerance = 1e-6
    penalty = 1e-3
    formulation = kinematic
    system = Constraint
  [../]
[]

[Materials]

  [./elasticity_tensor_Li]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 7.8e-3
    poissons_ratio = 0.3
    block = '2 8'
  [../]
  [./elasticity_tensor_Llzo]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 170e-3
    poissons_ratio = 0.25
    block = '1 5 '
  [../]
  [./stress_llzo]
    type = ADComputeFiniteStrainElasticStress
    block = '1 5 '
  [../]

  [./stress_Li]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    # perform_finite_strain_rotations = true
    block = '2'
  [../]

  [./stress_Li2]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas2'
    # perform_finite_strain_rotations = true
    block = '8'
  [../]
  [./plas]
    type = ADIsoTropicHyperVisco
    # absolute_tolerance = 1e-6
    # relative_tolerance = 1e-06
    hardening_exponent = 2.0
    saturation_resistance = 2.0e-6
    initial_resistance = 0.95e-6
    hardening_modulus = 10.0e-6
    rate_exponent = 0.15
    reference_strain_rate = 0.05
    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 0.1
    # omega = 2
    # alpha1 = 0.0
    # alpha2 = 1.0
    # alpha3 = 0.0
    # concentration = li_metal_conc
    cref = 0.0
    block = 2
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    enable = true
  [../]

  [./plas2]
    type = ADIsoTropicHyperVisco
    # relative_tolerance = 1e-06
    # absolute_tolerance = 1e-6
    hardening_exponent = 2.0
    saturation_resistance = 2.0e-6
    initial_resistance = 0.95e-6
    hardening_modulus = 10.0e-6
    rate_exponent = 0.15
    reference_strain_rate = 0.05
    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 0.1
    block = 8
  [../]

[]

[BCs]
  [./Li_top_y]
    type = ADPresetBC
    boundary = '18'
    variable = uy
    value = 0.0
  [../]

  [./Li_top_x]
    type = ADPresetBC
    boundary = '18'
    variable = ux
    value = 0.0
  [../]

  [./left_right]
    type = ADPresetBC
    variable = ux
    boundary = '14 17 12 19'
    value = 0
  [../]
  [./bot_pressure]
    type = ADPressure
    boundary = '11'
    variable = uy
    component = 1
    function = pressure
    # use_displaced_mesh = true
    # extra_vector_tags = 'ref'
    # type = ADFunctionDirichletBC
    # boundary = '11'
    # variable = uy
    # function = '-1.0*(t-0.01)'
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]

  type = Transient

  #Preconditioned JFNK (default)
  solve_type = 'NEWTON'
  automatic_scaling = true
  compute_scaling_once = false
  # petsc_options = '-snes_ksp_ew'
  # petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor'
  # petsc_options_iname = '-pc_type -pc_hypre_type -mat_mffd_err'
  # petsc_options_value = 'hypre    boomeramg      1e-8'
  # line_search = 'none'
  # petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart'
  # petsc_options_value = 'asm    lu     20 101'
  # petsc_options_iname = -pc_type
  # petsc_options_value = lu
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor'
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1     NONZERO               1e-20               '
  # line_search = 'contact'
  # line_search = 'none'
  # petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap
  #                       -ksp_gmres_restart'
  # petsc_options_value = 'asm lu 20 101'
  l_max_its = 30
  nl_max_its = 25




  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-10
  # l_tol = 1e-05

  # l_max_its = 100


  start_time = 0.0
  dt = 1.0
  dtmax = 0.5
  dtmin = 1e-5
  end_time = 200.0
  # num_steps = 10
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.025
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 40
    # timestep_limiting_postprocessor = matl_ts_min
  [../]

[] # Executioner

[Outputs]
  # [./check]
  #   type = Checkpoint
  #   num_files = 2
  #   interval = 5
  #   sync_only = true
  #   sync_times = '1 2.0 5.0 10.0 15.0 20.0 25.0 50.0 100.0 200.0'
  # [../]
  [./out]
    type = Exodus
    execute_on =  'INITIAL TIMESTEP_END'
    elemental_as_nodal = true
    # execute_elemental_variables = true
    execute_elemental_on = 'TIMESTEP_END'
    sync_only = false
    sync_times = '1.0 5.0 10.0 20.0 30.0 40.0 50.0 100.0 200.0'
    # sync_times = '1 2.0 5.0 10.0 15.0 20.0 25.0 50.0 100.0 200.0'
    # output_material_properties = true
  [../]
  [./csv]
    type = CSV
  [../]
[] # Outputs

[Debug]
  show_var_residual_norms = true
  show_material_props = true
[]
