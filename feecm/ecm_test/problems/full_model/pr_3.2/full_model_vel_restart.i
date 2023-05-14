elem = QUAD4
order = FIRST
[Mesh]
  patch_size = 80
  patch_update_strategy = auto
  parallel_type = REPLICATED
  [./mesh]
    type = FileMeshGenerator
    file = data/cosine_symmetric2.msh
  [../]
  [./secondary_boundary_block]
    type = LowerDBlockFromSidesetGenerator
    input = mesh
    sidesets = 'blockMetal_bottom'
    new_block_name = 'secondary'
  [../]
  [./primary_boundary_block]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_boundary_block
    sidesets = 'blockCeramic_top'
    new_block_name = 'primary'
  [../]
[]
[GlobalParams]
  displacements = 'ux uy'
[]

[Variables]
  [./ux]
    block = 'blockCeramic blockMetal interLayer'
  [../]
  [./uy]
    block = 'blockCeramic blockMetal interLayer'
  [../]
  [./li_ion_V]
    block = 'blockCeramic blockMetal interLayer'
  [../]
  [./thermal_lm]
    block = 'secondary'
  [../]
  [./normal_lm]
    block = 'secondary'
  [../]
  [./li_metal_conc]
    block = 'interLayer blockMetal'
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

  [./li_metal_flux_x]
    order = CONSTANT
    family = MONOMIAL
    block = 'blockCeramic interLayer blockMetal'
  [../]
  [./li_metal_flux_y]
    order = CONSTANT
    family = MONOMIAL
    block = 'blockCeramic interLayer blockMetal'
  [../]
  [./li_metal_flux_z]
    order = CONSTANT
    family = MONOMIAL
    block = 'blockCeramic interLayer blockMetal'
  [../]
  [./bndliflux]
    order = CONSTANT
    family = MONOMIAL
    block = 'blockCeramic interLayer blockMetal'
  [../]
  [./li_ion_flux_x]
    order = FIRST
    family = MONOMIAL
    block = 'blockCeramic interLayer blockMetal'
  [../]
  [./li_ion_flux_y]
    order = FIRST
    family = MONOMIAL
    block = 'blockCeramic interLayer blockMetal'
  [../]
  [./li_ion_flux_z]
    order = FIRST
    family = MONOMIAL
    block = 'blockCeramic interLayer blockMetal'
  [../]

[]
[AuxKernels]
  [./bnd_li_flux]
    type = DiffusionFluxNormalToBoundaryAux
    variable = bndliflux
    boundary = 'blockMetal_bottom blockCeramic_top'
    diffusion_variable = li_ion_V
    diffusivity = thermal_conductivity
  [../]

  [./peeq]
    type = ADMaterialRealAux
    variable = plastic_strain
    property = effective_plastic_strain
    execute_on = timestep_end
    block = 'interLayer blockMetal'
  [../]
  [./strain_rate]
    type = ADMaterialRealAux
    variable = strain_rate
    property = plastic_strain_rate
    block = 'interLayer blockMetal'
  [../]

  [./li_metal_flux_x]
    type = AnisoTropicDiffusionFluxAux
    variable = li_metal_flux_x
    component = x
    diffusion_variable = li_metal_conc
    block = 'interLayer blockMetal'
  [../]

  [./li_metal_flux_y]
    type = AnisoTropicDiffusionFluxAux
    variable = li_metal_flux_y
    component = y
    diffusion_variable = li_metal_conc
    block = 'interLayer blockMetal'
  [../]
  [./li_metal_flux_z]
    type = AnisoTropicDiffusionFluxAux
    variable = li_metal_flux_z
    component = z
    diffusion_variable = li_metal_conc
    block = 'interLayer blockMetal'
  [../]
  [./li_ion_flux_x]
    type = ADDiffusionFluxAux
    variable = li_ion_flux_x
    component = x
    diffusion_variable = li_ion_V
    diffusivity = thermal_conductivity
    block = 'blockCeramic blockMetal interLayer'
  [../]

  [./li_ion_flux_y]
    type = ADDiffusionFluxAux
    variable = li_ion_flux_y
    component = y
    diffusion_variable = li_ion_V
    diffusivity = thermal_conductivity
    block = 'blockCeramic blockMetal interLayer'
  [../]
  [./li_ion_flux_z]
    type = ADDiffusionFluxAux
    variable = li_ion_flux_z
    component = z
    diffusion_variable = li_ion_V
    diffusivity = thermal_conductivity
    block = 'blockCeramic blockMetal interLayer'
  [../]

[]
[Functions]
  [./pressure]
    # type = ParsedFunction
    # value = '0.1e-6'
    type = PiecewiseLinear
    data_file = 'pressure_time1.csv'
    format = columns
  [../]
  [./flux]
    type = ParsedFunction
    value = 'if (t >= 50.0, 4e-13, 0.0)'
  [../]
  [./gapk]
    # type = ParsedFunction
    # value = ' 1e-9*(1.0/(1.0 + exp(2.0*3e4*(x-1e-4))))'
    type = PiecewiseLinear
    data_file = 'gap_cond.csv'
    format = columns
  [../]
  [./gapk1]
    type = PiecewiseLinear
    x = '-10.0 0 1e-5 1e-4 1 10 100'
    y = '1e-2 1e-2 1e-2 0 0 0 0'
  [../]
  [./k_function]
    type = ParsedFunction
    value = 'if (t <=20, 1e-9, (1.0/(1.0 + exp(2.0*2.5*(x-6.0))))*(1e-9-1e-11) + 1e-11)'

  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress hydrostatic_stress'
    use_automatic_differentiation = true
    block = 'blockCeramic interLayer blockMetal'
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
    secondary_variable = li_ion_V
    primary_boundary =  'blockCeramic_top'
    secondary_subdomain = 'secondary'
    secondary_boundary = 'blockMetal_bottom'
    primary_subdomain = 'primary'
    k_function = gapk
    k = 1e-8
    use_displaced_mesh = true
    compute_lm_residuals = true
  [../]

  [./normal_lm]
    type = NormalNodalLMMechanicalContact
    secondary = 'blockMetal_bottom'
    primary = 'blockCeramic_top'
    variable = normal_lm
    primary_variable = ux
    normal_smoothing_distance = 0.2
    disp_y = uy
    ncp_function_type = min
    tangential_tolerance = 0.2
    c = 1e-3
  [../]
  [./normal_x]
    type = NormalMortarMechanicalContact
    primary_boundary = 'blockCeramic_top'
    secondary_boundary = 'blockMetal_bottom'
    primary_subdomain = 'primary'
    secondary_subdomain = 'secondary'
    variable = normal_lm
    secondary_variable = ux
    component = x
    use_displaced_mesh = true
    compute_lm_residuals = false
    # extra_vector_tags = 'ref'
  [../]
  [./normal_y]
    type = NormalMortarMechanicalContact
    secondary_boundary = 'blockCeramic_top'
    primary_boundary = 'blockMetal_bottom'
    primary_subdomain = 'primary'
    secondary_subdomain = 'secondary'
    variable = normal_lm
    secondary_variable = uy
    component = y
    use_displaced_mesh = true
    compute_lm_residuals = false
    # extra_vector_tags = 'ref'
  [../]
[]

[Kernels]
  [./li_ion_V]
    type = ADHeatConduction
    block = 'blockCeramic interLayer blockMetal'
    variable = li_ion_V
    use_displaced_mesh = false
  [../]
  [./li_metal2]
    type = ADMatAnisoDiffusion
    block = 'interLayer blockMetal'
    variable = li_metal_conc
    diffusivity = diffusivity
    use_displaced_mesh = false
  [../]
  [./li_metal_dt]
    type = ADTimeDerivative
    block = 'interLayer blockMetal'
    variable = li_metal_conc
    use_displaced_mesh = false
  [../]

[]

[Materials]
  [./diffusivity_Li]
    type = ADDiffusionAlongPrincipalDirections
    diffusivity_vector = '1e10 0 0'
    block = 'interLayer blockMetal'
  [../]
  [./thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1.0e-8
    block = 'blockCeramic'
  [../]

  [./thermal_conductivity2]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1e-4
    block = 'interLayer blockMetal'
    use_displaced_mesh = true
  [../]

  [./elasticity_tensor_Li]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 7.8e-3
    poissons_ratio = 0.3
    block = 'interLayer blockMetal'
  [../]
  [./elasticity_tensor_Llzo]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 170e-3
    poissons_ratio = 0.25
    block = 'blockCeramic'
  [../]
  [./stress_llzo]
    type = ADComputeFiniteStrainElasticStress
    block = 'blockCeramic'
  [../]

  [./stress_Li]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    # perform_finite_strain_rotations = true
    block = 'interLayer'
  [../]

  [./stress_Li2]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas2'
    # perform_finite_strain_rotations = true
    block = 'blockMetal'
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
    temperature = 333

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 0.1
    omega = 2
    alpha = '1 0 0'
    concentration = li_metal_conc
    cref = 0
    block = 'interLayer'
    intBnd = 'blockMetal_bottom'
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    enable = true
  [../]

  [./plas2]
    type = ADIsoTropicHyperViscoSwellingCreep
    hardening_exponent = 2.0
    saturation_resistance = 2.0e-6
    initial_resistance = 0.95e-6
    hardening_modulus = 10.0e-6
    rate_exponent = 0.15
    activation_energy = 37
    gas_constant = 8.314462681e-3 # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = 333

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 0.1
    omega = 2e-6
    alpha = '1 0 0'
    intBnd = 'blockMetal_top'
    concentration = li_metal_conc
    cref = 0
    block = 'blockMetal'
  [../]
[]

[BCs]
  [./Li_top_y]
    type = ADDirichletBC
    variable = uy
    boundary = 'blockMetal_top'
    value = 0.0
  [../]

  [./Li_top_x]
    type = ADDirichletBC
    boundary = 'blockMetal_top'
    variable = ux
    value = 0.0
  [../]


  [./left_right]
    type = ADDirichletBC
    variable = ux
    boundary = 'blockCeramic_left blockCeramic_right blockMetal_left blockMetal_right'
    value = 0
  [../]
  [./li_ion_bottom]
    type = FunctionNeumannBC
    boundary = 'blockCeramic_bottom'
    variable = li_ion_V
    function = flux
    extra_vector_tags = 'ref'
  [../]

  [./li_ion_top]
    type = ADDirichletBC
    boundary = 'blockMetal_top'
    variable = li_ion_V
    value = 0.0
  [../]

  [./top_pressure]
    type = PresetVelocity
    boundary = 'blockCeramic_bottom'
    variable = uy
    # component = 1
    # function = pressure
    displacements = 'ux uy'
    # use_displaced_mesh = true
    velocity = 1e-6
    function = 'if (t <= 50.0, 1,0.0)'
  [../]



  [./li_metal_flux]
    type = ScaledCoupledVarNeumannBC
    boundary = 'blockMetal_bottom'
    variable = li_metal_conc
    v = bndliflux
    scale = -1.036426e9
    extra_vector_tags = 'ref'
    # value = 1.0
  [../]

[]

[Postprocessors]
  [./bottom_current]
    type = ADSideFluxIntegral
    variable = li_ion_V
    boundary = 'blockCeramic_top'
    diffusivity = thermal_conductivity
  [../]
  [./top_current]
    type = ADSideFluxIntegral
    variable = li_ion_V
    boundary = 'blockMetal_bottom'
    diffusivity = thermal_conductivity
  [../]
  [./lagrange]
    type = ElementIntegralVariablePostprocessor
    variable = thermal_lm
    block = 'secondary'
  [../]
  # [./auxcurrent]
  #   type = SideIntegralVariablePostprocessor
  #   variable = bndliflux
  #   boundary = 15
  # [../]
  [./over_potential]
    type = SideAverageValue
    variable = li_ion_V
    boundary = 'blockCeramic_bottom'
  [../]

  [./top_pressure]
    type = SideAverageValue
    variable = stress_yy
    boundary = 'blockMetal_top'
  [../]
  [./bot_pressure]
    type = SideAverageValue
    variable = stress_yy
    boundary = 'blockCeramic_bottom'
  [../]
  [./matl_ts_min]
    type = MaterialTimeStepPostprocessor
    block = 'interLayer blockMetal'
  [../]
  [./contact_size]
    type = ContactDOFSetSize
    variable = normal_lm
    subdomain = 'secondary'
    tolerance = 1e-12
  [../]

  # [./limetal_flux]
  #   type = SideFluxIntegral
  #   variable = li_metal_conc
  #   boundary = 15
  #   diffusivity = diffusivity_Li
  # [../]
[]

# [VectorPostProcessor]
#   [./flux_interface]
#     type =
#   [../]
# []

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
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor'
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1     NONZERO               1e-20               '
  line_search = contact
  l_max_its = 50
  nl_max_its = 25
  nl_abs_tol = 1e-16
  nl_rel_tol = 1e-9
  start_time = 0.0
  dt = 1.0
  dtmax = 0.5
  dtmin = 1e-5
  end_time = 500.0
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
  [./out]
    type = Exodus
    sync_times = '9.9 10.0 20.0 25.0 30.0 50.0 80.0 110.0 150.0 180.0 210.0 240.0 270.0 300.0 330.0 360.0 500.0'
    sync_only = false
    # interval = 2
  [../]
  [./csv]
    type = CSV
  [../]
  [./checkpoint]
    type = Checkpoint
    num_files = 4
    sync_times = '9.9 10.0 20.0 40.0 50.0 70.0 100.0 200.0 250.0 300.0 400.0 500.0'
    sync_only = true
  [../]

[] # Outputs

[Debug]
  show_var_residual_norms = true
  show_material_props = true
[]
