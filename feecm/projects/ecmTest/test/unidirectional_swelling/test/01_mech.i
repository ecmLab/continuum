# This is a model with 2 flat blocks
# block 1 - Ceramic Electrolyte (top)
# block 2 - Li_metal with Anand model swelling
# Li_ion_conc -- modeled as a steady state thermal problem
# All units in um
elem = QUAD4
order = FIRST
name = 'mech'

[Mesh]
  patch_size = 80
  patch_update_strategy = auto
  displacements = 'disp_x disp_y'
#  parallel_type = REPLICATED
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

[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'disp_x disp_y'
  acceptable_iterations = 2
[]

[Variables]
  [./disp_x]
    block = '1 2 5 8'
    order = ${order}
    # scaling = 1e12
  [../]
  [./disp_y]
    block = '1 2 5 8'
    order = ${order}
    # scaling = 1e12
  [../]

[]

[Functions]
  [./pressure]
    type = ParsedFunction
    value = 'if (t <= 10.0, 1.0e-5*t + 0.4e-5, 1.4e-5)'
  [../]
  [./fl_x]
    type = ParsedFunction
    value = 'if (t > 10, if (t>= 20, -4e-12, -4e-12*(t-10)/10.0), 5e-10)'
  [../]
  [./gapk]
    type = PiecewiseLinear
    x = '-10.0 0 1e-5 1e-4 1 10 100'
    y = '1e-9 1e-9 1e-9 0 0 0 0'
  [../]
  [./gapk1]
    type = PiecewiseLinear
    x = '-10.0 0 1e-5 1e-4 1 10 100'
    y = '1e-2 1e-2 1e-2 0 0 0 0'
  [../]
  [./k_function]
    type = ParsedFunction
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
    displacements = 'disp_x disp_y'
    extra_vector_tags = 'ref'
  [../]
[]

[Contact]
  [./mech]
    primary = 16
    secondary = 15
    displacements = 'disp_x disp_y'
    capture_tolerance = 1e-6
    penalty = 1e-3
    mesh = cmg
    formulation = mortar
#    formulation = kinematic
  [../]
[]

[BCs]
  [./SE_top_y]
    type = ADPresetBC
    boundary = '18'
    variable = disp_y
    value = 0.0
  [../]

  [./left_right]
    type = ADPresetBC
    variable = disp_x
    boundary = '12 14 17 18 19'
    value = 0
  [../]
  [./bot_pressure]
    type = ADPressure
    boundary = '11'
    variable = disp_y
    component = 1
    function = pressure
  [../]
[]

[Materials]
  [./elasticity_tensor_Li]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 7.8e-3
    poissons_ratio = 0.38
    block = '2 8'
  [../]
  [./elasticity_tensor_Llzo]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 170e-3
    poissons_ratio = 0.3
    block = '1 5'
  [../]
  [./stress_llzo]
    type = ADComputeFiniteStrainElasticStress
    block = '1 5'
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
    type = ADIsoTropicHyperViscoSwellingCreep
    # absolute_tolerance = 1e-6
    # relative_tolerance = 1e-06
    hardening_exponent = 2.0
    saturation_resistance = 2.0e-6
    initial_resistance = 0.95e-6
    hardening_modulus = 10.0e-6
    rate_exponent = 0.15
    # reference_strain_rate = 0.05
    activation_energy = 37
    gas_constant = 8.314462681e-3 # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = 298

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 0.1
    omega = 2
    alpha1 = 0.0
    alpha2 = 1.0
    alpha3 = 0.0
    concentration = 0
    cref = 0.0
    block = 2
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    enable = true

  [../]

  [./plas2]
    type = ADIsoTropicHyperViscoCreep
    # relative_tolerance = 1e-06
    # absolute_tolerance = 1e-6
    hardening_exponent = 2.0
    saturation_resistance = 2.0e-6
    initial_resistance = 0.95e-6
    hardening_modulus = 10.0e-6
    rate_exponent = 0.15
    # reference_strain_rate = 0.05
    activation_energy = 37
    gas_constant = 8.314462681e-3 # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = 298

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 0.1
    block = 8
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
#  compute_scaling_once = false
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor'
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1     NONZERO               1e-20               '

  line_search = 'none'
  nl_abs_tol = 1e-10
  l_max_its = 500
  nl_max_its = 30

  dt = 0.1
  dtmax = 2.0
  dtmin = 1e-5
  end_time = 20.0
  # num_steps = 10
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 40
    # timestep_limiting_postprocessor = matl_ts_min
  [../]

[] # Executioner

[Outputs]
  exodus = true
  csv = true
  file_base = rst/${name}
  elemental_as_nodal = true
  [out]
    type = Checkpoint
    interval = 5
    num_files = 2
  []
[]

[Debug]
  show_var_residual_norms = true
  show_material_props = true
[]
