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
    ymax = 50
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
  [./slave_boundary_block]
    type = LowerDBlockFromSidesetGenerator
    input = cmg
    sidesets = '16'
    new_block_id = '4'
  [../]
  [./master_boundary_block]
    type = LowerDBlockFromSidesetGenerator
    input = slave_boundary_block
    sidesets = '15'
    new_block_id = '3'
  [../]

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
  [./li_ion_V]
    block = '1 2 5 8'
    order = ${order}
    # scaling = 1e8
  [../]
  [./thermal_lm]
    block = '3'
    order = ${order}
  [../]
  # [./thermal_lm2]
  #   block = '3'
  #   order = ${order}
  # [../]

  [./normal_lm]
    block = 3
    order = ${order}
  [../]

  [./li_metal_conc]
    block = '2 8'
    initial_condition = 0
    order = ${order}
    # scaling = 1e3
  [../]

[]

[AuxVariables]

  [./plastic_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_rate]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./swelling]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./li_ion_flux_x]
    order = FIRST
    family = MONOMIAL
    block = '1 2 5 8'
  [../]
  [./li_ion_flux_y]
    order = FIRST
    family = MONOMIAL
    block = '1 2 5 8'
  [../]
  [./li_ion_flux_z]
    order = FIRST
    family = MONOMIAL
    block = '1 2 5 8'
  [../]
  [./bndliflux]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 5 8'
  [../]

[]
[AuxKernels]
  [./bnd_li_flux]
    type = DiffusionFluxNormalToBoundaryAux
    variable = bndliflux
    boundary = '15 16'
    diffusion_variable = li_ion_V
    diffusivity = thermal_conductivity
  [../]

[Outputs]
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
    file_base = rst/${name}
  [../]
  [./csv]
    type = CSV
  [../]
[] # Outputs

  [./peeq]
    type = ADMaterialRealAux
    variable = plastic_strain
    property = effective_plastic_strain
    execute_on = timestep_end
    block = '2 8'
  [../]
  [./strain_rate]
    type = ADMaterialRealAux
    variable = strain_rate
    property = plastic_strain_rate
    block = '2 8'
  [../]
  [./swelling]
    type = ADMaterialRealAux
    variable = swelling
    property = swelling_strain
    block = '2'
  [../]
  [./li_ion_flux_x]
    type = ADDiffusionFluxAux
    variable = li_ion_flux_x
    component = x
    diffusion_variable = li_ion_V
    diffusivity = thermal_conductivity
    block = '1 2 5  8'
  [../]

  [./li_ion_flux_y]
    type = ADDiffusionFluxAux
    variable = li_ion_flux_y
    component = y
    diffusion_variable = li_ion_V
    diffusivity = thermal_conductivity
    block = '1 2 5 8'
  [../]
  [./li_ion_flux_z]
    type = ADDiffusionFluxAux
    variable = li_ion_flux_z
    component = z
    diffusion_variable = li_ion_V
    diffusivity = thermal_conductivity
    block = '1 2 5 8'
  [../]

[]

[Functions]
  [./pressure]
    # type = ParsedFunction
    # value = '1.2e-6'
    type = ParsedFunction
    value = 'if (t <= 60.0, 0.05e-6*t + 0.2e-6, 2.0e-6)'
    # type = PiecewiseLinear
    # data_file = 'pressure_time.csv'
    # format = columns
    # direction = left
  [../]
  [./flux]
    type = ParsedFunction
    value = 'if (t > 60, if (t>= 70, -4e-12, -4e-12*(t-60)/9.0), 8e-11)'
    # value = '-2e-12/(1.0 + exp(-2.0*0.6*(t-5.0)))'
  [../]
  [./gapk]
    # type = ParsedFunction
    # value = ' 1e-9*(1.0/(1.0 + exp(2.0*3e4*(x-1e-4))))'
    type = PiecewiseLinear
    x = '-10.0 0 1e-6 1e-5 1 10 100'
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
    value = 'if (t <=70, 1e-8, (1.0/(1.0 + exp(2.0*2.5*(x-3.0))))*(1e-9-1e-11) + 1e-11)'

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
  # solution_variables = 'ux uy normal_lm thermal_lm li_ion_V thermal_lm2 li_metal_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'ux uy normal_lm; li_ion_V li_metal_conc thermal_lm'
  acceptable_iterations = 2
[]

[Constraints]
  [./thermal_constraint]
    type = GapDisplacementConductanceConstraint
    variable = thermal_lm
    slave_variable = li_ion_V
    slave_boundary = 15
    slave_subdomain = 3
    master_boundary = 16
    master_subdomain = 4
    k_function = gapk
    k = 1e-8
    use_displaced_mesh = true
    compute_lm_residuals = true
    # extra_vector_tags = 'ref'
  [../]

  # [./thermal_constraint1]
  #   type = GapDisplacementConductanceConstraint
  #   variable = thermal_lm2
  #   slave_variable = li_metal_conc
  #   slave_boundary = 15
  #   slave_subdomain = 3
  #   master_boundary = 16
  #   master_subdomain = 4
  #   k_function = gapk1
  #   k = 1e-8
  #   use_displaced_mesh = true
  #   compute_lm_residuals = true
  #   # extra_vector_tags = 'ref'
  # [../]

  [./normal_lm]
    type = NormalNodalLMMechanicalContact
    slave = 15
    master = 16
    variable = normal_lm
    master_variable = ux
    normal_smoothing_distance = 0.2
    disp_y = uy
    ncp_function_type = min
    tangential_tolerance = 0.2
    c = 1e-3
    # c = 1e-3
    # extra_vector_tags = 'ref'
    # c = 1e3
    # ncp_function_type = min
    # order = FIRST
    # tangential_tolerance = 0.1
  [../]
  [./normal_x]
    type = NormalMortarMechanicalContact
    master_boundary = 16
    slave_boundary = 15
    master_subdomain = 4
    slave_subdomain = 3
    variable = normal_lm
    slave_variable = ux
    component = x
    use_displaced_mesh = true
    compute_lm_residuals = false
    # extra_vector_tags = 'ref'
  [../]
  [./normal_y]
    type = NormalMortarMechanicalContact
    master_boundary = 16
    slave_boundary = 15
    master_subdomain = 4
    slave_subdomain = 3
    variable = normal_lm
    slave_variable = uy
    component = y
    use_displaced_mesh = true
    compute_lm_residuals = false
    # extra_vector_tags = 'ref'
  [../]

[]
[Kernels]
  [./li_ion_V]
    type = ADHeatConduction
    block = '1 2 5 8'
    variable = li_ion_V
    use_displaced_mesh = false
  [../]
  # [./li_metal1]
  #   type = ADMatDiffusion
  #   block = '1 5'
  #   variable = li_metal_conc
  #   diffusivity = diffusivity_Li1
  #   # extra_vector_tags = 'ref2'
  #   use_displaced_mesh = false
  # [../]
  [./li_metal2]
    type = ADMatAnisoDiffusion
    block = '2 8'
    variable = li_metal_conc
    diffusivity = diffusivity_Li
    # extra_vector_tags = 'ref'
    use_displaced_mesh = false
  [../]
  [./li_metal_dt]
    type = ADTimeDerivative
    block = '2 8'
    # extra_vector_tags = 'ref'
    # extra_vector_tags = 'ref2'
    variable = li_metal_conc
    use_displaced_mesh = false
  [../]

[../]


[Materials]
  # [./diffusivity_Li]
  #   type = GenericConstantMaterial
  #   block = '1 5'
  #   prop_names = 'diffusivity_Li1 specific_heat1'
  #   prop_values = '1e-3 1'
  # [../]
  # [./diffusivity_Li1]
  #   type = GenericConstantMaterial
  #   block = '2 8'
  #   prop_names = 'diffusivity_Li specific_heat1'
  #   prop_values = '1e2 1'
  # [../]

  [./diffusivity_Li1]
    # type = GenericConstantMaterial
    # prop_names = 'diffusivity_Li specific_heat1'
    # prop_values = '1e-1 1'
    block = '2 8'
    type = ADConstantAnisotropicMobility
    tensor = '0 0 0
              0 1e15 0
              0 0 0'
    M_name = 'diffusivity_Li'
  [../]

  [./thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1.0e-8
    block = '1'
  [../]

  [./thermal_function]
    type = ADGenericFunctionMaterial
    prop_names = 'thermal_conductivity'
    prop_values = k_function
    block = '5'
  [../]

  [./thermal_conductivity4]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1e-4
    block = '2 8'
    use_displaced_mesh = true
  [../]

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
    concentration = li_metal_conc
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
  [./li_ion_top]
    type = FunctionNeumannBC
    boundary = '18'
    variable = li_ion_V
    function = flux
    extra_vector_tags = 'ref'
  [../]

  [./li_ion_bot]
    type = ADPresetBC
    boundary = '11'
    variable = li_ion_V
    value = 0.0
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
  [./li_metal_flux]
    type = ScaledCoupledVarNeumannBC
    boundary = '15'
    variable = li_metal_conc
    v = bndliflux
    scale = -1.036426e9
    extra_vector_tags = 'ref'
    # value = 1.0
  [../]
  # [./li_con_top]
  #   type = ADPresetBC
  #   boundary = 18
  #   variable = li_metal_conc
  #   value = 0
  # [../]

[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./bottom_current]
    type = ADSideFluxIntegral
    variable = li_ion_V
    boundary = 16
    diffusivity = thermal_conductivity
  [../]
  [./top_current]
    type = ADSideFluxIntegral
    variable = li_ion_V
    boundary = 15
    diffusivity = thermal_conductivity
  [../]
  [./lagrange]
    type = ElementIntegralVariablePostprocessor
    variable = thermal_lm
    block = 3
  [../]
  # [./auxcurrent]
  #   type = SideIntegralVariablePostprocessor
  #   variable = bndliflux
  #   boundary = 15
  # [../]
  [./over_potential]
    type = SideAverageValue
    variable = li_ion_V
    boundary = 18
  [../]

  [./ext_pressure]
    type = SideAverageValue
    variable = stress_yy
    boundary = '11'
  [../]
  [./matl_ts_min]
    type = MaterialTimeStepPostprocessor
    block = '2 8'
  [../]

  # [./limetal_flux]
  #   type = SideFluxIntegral
  #   variable = li_metal_conc
  #   boundary = 15
  #   diffusivity = diffusivity_Li
  # [../]
[]
# [VectorPostprocessors]
#   [./boundary_flux]
#     type = SideValueSampler
#     boundary = 15
#     variable = bndliflux
#     sort_by = x
#     contains_complete_history = true
#   [../]
# []

[Executioner]

  type = Transient

  #Preconditioned JFNK (default)
  solve_type = 'NEWTON'
  automatic_scaling = true
  compute_scaling_once = false
  # petsc_options = '-snes_ksp_ew'
  # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  # petsc_options_value = 'hypre    boomeramg      101'
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
  l_max_its = 50
  nl_max_its = 25




  nl_abs_tol = 3e-15
  nl_rel_tol = 1e-10
  # l_tol = 1e-05

  # l_max_its = 100


  start_time = 0.0
  dt = 1.0
  dtmax = 2.0
  dtmin = 1e-5
  end_time = 500.0
  # num_steps = 10
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.025
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 40
    timestep_limiting_postprocessor = matl_ts_min
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
    sync_times = '1.0 5.0 10.0 20.0 30.0 40.0 50.0 100.0 200.0 250.0 300.0 400.0 500.0'
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

# [Controls]
#   [./kfunc]
#     type = TimePeriod
#     start_time = 0
#     end_time = 1
#     enable_objects = 'constkfunc'
#     disable_objects = 'statkfunc'
#   [../]
# []
