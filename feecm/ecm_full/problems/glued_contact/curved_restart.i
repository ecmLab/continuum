[Mesh]
  patch_update_strategy = iteration
  patch_size = 200
  parallel_type = REPLICATED
  file = check/curved_cp/0900_mesh.cpr
  [./mesh]
    type = FileMeshGenerator
    file = data/initialDefect2.msh
  [../]
  # [./secondary_boundary_block]
  #   type = LowerDBlockFromSidesetGenerator
  #   input = mesh
  #   sidesets = 'blockMetal_bottom'
  #   new_block_name = 'mech_contact_primary_subdomain'
  # [../]
  # [./primary_boundary_block]
  #   type = LowerDBlockFromSidesetGenerator
  #   input = secondary_boundary_block
  #   sidesets = 'blockCeramic_top'
  #   new_block_name = 'mech_contact_secondary_subdomain'
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
  [./li_ion_V]
    block = 'blockCeramic blockMetal interLayer'
  [../]
  [./thermal_lm]
    block = 'mech_contact_secondary_subdomain'
  [../]
  [./li_metal_conc2]
    block = 'interLayer'
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
  [./bndliflux]
    order = CONSTANT
    family = MONOMIAL
    block = 'blockCeramic interLayer blockMetal'
  [../]

  [./li_metal_conc]
    block = 'interLayer blockMetal'
    # initial_condition = 0
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
  [./swelling]
    order = CONSTANT
    family = MONOMIAL
    block = 'interLayer blockMetal'
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
  [./li_metal_conc]
    type = ConstantAux
    value = 0
    variable = li_metal_conc
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
  [./swelling_vol_change]
    type = ADMaterialRealAux
    variable = swelling
    block = 'interLayer blockMetal'
    property = swelling_vol_change
  [../]


[]

[Functions]
  [./pressure]
    type = ParsedFunction
    value = '1.0e-5'
  [../]
  [./uy]
    type = ParsedFunction
    value = '0.01*t'
  [../]
  [./vel2]
    type = ParsedFunction
    value = '1e-9*t'
  [../]
  [./flux]
    type = ParsedFunction
    value = 'if (t <= 1, if (t > 0.001, 1e-13 + 9e-13*t,0.0), 1e-12)'
  [../]
  [./gapk]
    # type = ParsedFunction
    # value = ' 1e-9*(1.0/(1.0 + exp(2.0*3e4*(x-1e-4))))'
    type = PiecewiseLinear
    data_file = data/gap_cond.csv
    format = columns
  [../]
  [./gapk1]
    type = PiecewiseLinear
    x = '-10.0 0 1e-5 1e-4 1 10 100'
    y = '1e-2 1e-2 1e-2 0 0 0 0'
  [../]

[]

[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress'
    use_automatic_differentiation = true
    extra_vector_tags = 'ref'
    block = 'blockMetal blockCeramic interLayer'
    # use_finite_deform_jacobian = true
  [../]
[]
[Problem]
  type = ReferenceResidualProblem
  # solution_variables = 'ux uy normal_lm thermal_lm li_ion_V thermal_lm2 blockMetal_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'ux uy mech_contact_normal_lm; li_ion_V thermal_lm li_metal_conc2'
  acceptable_iterations = 2
  restart_file_base = check/curved_cp/0900
[]

[Contact]
  [./mech_contact]
    disp_x = ux
    disp_y = uy
    secondary  =  'blockCeramic_top'
    primary =  'blockMetal_bottom'
    # penalty = 1e-1
    # formulation = penalty
    normal_smoothing_distance = 0.5
    tangential_tolerance = 0.25
    formulation = mortar
    mesh = mesh
    use_dual = true
  [../]
[]
[Constraints]
  [./thermal_constraint]
    type = GapDisplacementConductanceConstraint
    variable = thermal_lm
    secondary_variable = li_ion_V
    secondary_boundary =  'blockCeramic_top'
    secondary_subdomain = 'mech_contact_secondary_subdomain'
    primary_boundary = 'blockMetal_bottom'
    primary_subdomain = 'mech_contact_primary_subdomain'
    k_function = gapk
    k = 1e-8
    pressure = mech_contact_normal_lm
    use_displaced_mesh = true
    compute_lm_residuals = true
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
    block = 'interLayer'
    variable = li_metal_conc2
    diffusivity = diffusivity
    use_displaced_mesh = false
  [../]
  [./li_metal_dt]
    type = ADTimeDerivative
    block = 'interLayer'
    variable = li_metal_conc2
    use_displaced_mesh = false
  [../]

[]

[Materials]
  [./diffusivity_Li]
    type = ADDiffusionAlongPrincipalDirections
    diffusivity_vector = '1e2 0 0'
    block = 'interLayer'
  [../]

 # [./diffusivity_Li2]
 #   type = ADDiffusionAlongPrincipalDirections
 #   diffusivity_vector = '1 0 0'
 #   block = 'blockMetal'
 # [../]

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
    youngs_modulus = 7.81e-3
    poissons_ratio = 0.38
    block = 'blockMetal interLayer'
  [../]
  [./elasticity_tensor_Llzo]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 170e-3
    poissons_ratio = 0.3
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

  [./plas]
    type = ADIsoTropicHyperViscoSwellingCreep
    # absolute_tolerance = 1e-6
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
    alpha = '1 0 0'
    intBnd = 'blockCeramic_top'
    concentration = li_metal_conc2
    cref = 0
    omega = 2
    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 1.0
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    block = 'interLayer'
  [../]

  [./stress_Li2]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas2'
    # perform_finite_strain_rotations = true
    block = 'blockMetal'
  [../]

  [./plas2]
    type = ADIsoTropicHyperViscoSwellingCreep
    # absolute_tolerance = 1e-6
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
    alpha = '1 0 0'
    intBnd = 'blockMetal_top'
    concentration = li_metal_conc
    cref = 0
    omega = 0.0
    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 1.0
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    block = 'blockMetal'
  [../]

[]

[BCs]
  [./pressure]
    type = PresetVelocity
    variable = uy
   component = 1
    function = vel2
    boundary = 'blockCeramic_bottom'
  [../]

  [./left_right]
    type = ADDirichletBC
    preset =  true
    variable = ux
    boundary = 'blockCeramic_left blockCeramic_bottom blockCeramic_right blockMetal_right blockMetal_left blockMetal_top'
    value = 0
  [../]
  [./bottom]
    type = ADDirichletBC
    preset = true
    variable = uy
    boundary = 'blockMetal_top'
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

  [./li_metal_flux]
    type = ScaledCoupledVarNeumannBC
    boundary = 'blockMetal_bottom'
    variable = li_metal_conc2
    v = bndliflux
    scale = -1.036426e9
    extra_vector_tags = 'ref'
    # value = 1.0
  [../]

  # [./bottom_flux]
  #   type = ADNeumannBC
  #   # preset = true
  #   variable = li_metal_conc
  #   boundary = 'blockMetal_top'
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
  [./matl_ts_min]
    type = MaterialTimeStepPostprocessor
    block = 'interLayer blockMetal'
  [../]
  [./bot_stress]
    type = SideAverageValue
    variable = stress_yy
    boundary = 'blockCeramic_bottom'
  [../]
  [./bottom_interface_current]
    type = ADSideFluxIntegral
    variable = li_ion_V
    boundary = 'blockCeramic_top'
    diffusivity = thermal_conductivity
  [../]
  [./top_interface_current]
    type = ADSideFluxIntegral
    variable = li_ion_V
    boundary = 'blockMetal_bottom'
    diffusivity = thermal_conductivity
  [../]
  [./lagrange]
    type = ElementIntegralVariablePostprocessor
    variable = thermal_lm
    block = 'mech_contact_secondary_subdomain'
  [../]

  [./over_potential]
    type = SideAverageValue
    variable = li_ion_V
    boundary = 'blockCeramic_bottom'
  [../]

  [./ext_pressure]
    type = SideAverageValue
    variable = stress_yy
    boundary = 'blockMetal_top'
  [../]
  [./contact_size]
    type = ContactDOFSetSize
    variable = mech_contact_normal_lm
    subdomain = 'mech_contact_secondary_subdomain'
    tolerance = 1e-15
  [../]

[]

[VectorPostprocessors]
  [./bound_flux_ceramic]
    type = SideValueSampler
    variable = 'li_ion_flux_x li_ion_flux_y'
    boundary = 'blockCeramic_top'
    sort_by = x
  [../]
  [./bound_flux_interlayer]
    type = SideValueSampler
    variable = 'li_ion_flux_x li_ion_flux_y'
    boundary = 'blockMetal_bottom'
    sort_by = x
  [../]

[]

[Executioner]

  type = Transient

  #Preconditioned JFNK (default)
  solve_type = 'NEWTON'
  automatic_scaling = true
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor -snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_force_iteration -pc_factor_shift_type -pc_factor_shift_amount -mat_mffd_err'
  petsc_options_value = 'lu superlu_dist 1     NONZERO               1e-15               1e-5'



  line_search = 'none'


  nl_abs_tol = 3e-15
  nl_rel_tol = 1e-10

  l_max_its = 100
  nl_max_its = 35

  start_time = 0.0
  dt = 0.001
  dtmax = 0.05
  dtmin = 1e-5
  end_time = 200.0
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.001
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 30
    timestep_limiting_postprocessor = matl_ts_min
  [../]

[] # Executioner

[Outputs]
  [./out]
    type = Exodus
    # elemental_as_nodal = true
    # sync_times = '0.1 0.2 0.5 0.6 0.7 0.8 0.9 1.0 1.5 2.0'
    sync_only = false
    file_base = rst/curved_restart
  [../]
  [./csv]
    type = CSV
    file_base = rst/curved_restart
  [../]
  [./check]
    type = Checkpoint
    num_files = 10 
    start_time = 0.0
    sync_only = false
    sync_times = '0.005 0.2 0.5 1.0 2.0 4.0 5.0 8.0 9.0 10.0 20.0 25.0 30.0 40.0 100.0 150.0 200.0'
    file_base = check/curved_restart
    interval = 50
  [../]
[] # Outputs
