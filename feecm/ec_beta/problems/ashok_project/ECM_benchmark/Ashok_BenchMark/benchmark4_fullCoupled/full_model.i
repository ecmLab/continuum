# elem = QUAD4
# order = FIRST
[Mesh]
  patch_size = 80
  patch_update_strategy = auto
  parallel_type = REPLICATED
  # file = check/full_model_cp/0002_mesh.cpr
  [./mesh]
    type = FileMeshGenerator
    file = data/cosine_conformal2.msh
  [../]
  # [./primary_block]
  #   type = LowerDBlockFromSidesetGenerator
  #   input = mesh
  #   sidesets = 'blockMetal_bottom'
  #   new_block_name = 'mech_contact_primary_subdomain'
  # [../]
  # [./secondary_block]
  #   type = LowerDBlockFromSidesetGenerator
  #   input = primary_block
  #   sidesets = 'blockCeramic_top'
  #   new_block_name = 'mech_contact_secondary_subdomain'
  # [../]

[]
[Problem]
  coord_type = 'RZ'
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
    block = 'mech_contact_secondary_subdomain'
  [../]

  [./li_metal_conc]
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
    boundary = 'blockCeramic_top blockMetal_bottom'
    diffusion_variable = li_ion_V
    diffusivity = thermal_conductivity
    check_boundary_restricted = false
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
    block = 'interLayer'
  [../]

  [./li_metal_flux_y]
    type = AnisoTropicDiffusionFluxAux
    variable = li_metal_flux_y
    component = y
    diffusion_variable = li_metal_conc
    block = 'interLayer'
  [../]
  [./li_metal_flux_z]
    type = AnisoTropicDiffusionFluxAux
    variable = li_metal_flux_z
    component = z
    diffusion_variable = li_metal_conc
    block = 'interLayer'
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
    type = ParsedFunction
    value = 'if (t<=1,0.2 + 1.8*t,2.0)'
    # type = PiecewiseLinear
    # data_file = 'pressure_time1.csv'
    # format = columns
  [../]
  [./flux]
    type = ParsedFunction
    value = 'if (t > 1.0, if (t < 6.0, 0.2e-4*(t-1), 1e-3), 0.0)'
  [../]
  [./gapk]
    type = PiecewiseLinear
    x = '0 1e-15 1e-14 1    1e6   1e9'
    y = '0  0    1e-3  1e-2 1e-3    1e-3 '
  [../]
  [./gapk1]
    type = PiecewiseLinear
    x = '-5   -1e-3 -1e-6 -1e-15 1e-12 1e-11 1'
    y = '1e-3  1e-3  1e-3   1e-3  1e-3  0    0'
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
  group_variables = 'ux uy mech_contact_normal_lm; li_ion_V li_metal_conc thermal_lm'
  acceptable_iterations = 2
  # restart_file_base = check/full_model_cp/0002
[]

[Contact]
  [./mech_contact]
    # disp_x = ux
    # disp_y = uy
    mesh = mesh
    secondary = 'blockCeramic_top'
    primary = 'blockMetal_bottom'
    normal_smoothing_distance = 0.5
    tangential_tolerance = 0.25
    formulation = mortar
  [../]
[]

[Constraints]
  [./thermal_constraint]
    type = PressureLMConductanceConstraint
    contact_pressure = mech_contact_normal_lm
    variable = thermal_lm
    secondary_variable = li_ion_V
    secondary_boundary =  'blockCeramic_top'
    secondary_subdomain = 'mech_contact_secondary_subdomain'
    primary_boundary = 'blockMetal_bottom'
    primary_subdomain = 'mech_contact_primary_subdomain'
    k = 1e-3
    use_displaced_mesh = true
    # compute_lm_residuals = true
    # compute_primal_residuals = true
  [../]

  [./therm2]
    type = ScaledBCConstraint
    variable = thermal_lm
    primary_variable = li_metal_conc
    secondary_variable = li_ion_V
    secondary_boundary =  'blockCeramic_top'
    secondary_subdomain = 'mech_contact_secondary_subdomain'
    primary_boundary = 'blockMetal_bottom'
    primary_subdomain = 'mech_contact_primary_subdomain'
    scale = -1.036426e-2
    use_displaced_mesh = true
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
    variable = li_metal_conc
    diffusivity = diffusivity
    use_displaced_mesh = false
  [../]
  [./li_metal_dt]
    type = ADTimeDerivative
    block = 'interLayer'
    variable = li_metal_conc
    use_displaced_mesh = false
  [../]

[]

[Materials]
  [./diffusivity_Li]
    type = ADDiffusionAlongPrincipalDirectionsMaterial   #ADDiffusionAlongPrincipalDirections
    diffusivity_vector = '1e10 0 0'
    block = 'interLayer'
  [../]
  [./thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1.0e-2
    block = 'blockCeramic'
  [../]

  [./thermal_conductivity2]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1e2
    block = 'interLayer blockMetal'
    use_displaced_mesh = true
  [../]

  [./elasticity_tensor_Li]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 7.8e3
    poissons_ratio = 0.3
    block = 'interLayer blockMetal'
  [../]
  [./elasticity_tensor_Llzo]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 170e3
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
    absolute_tolerance = 1e-8
    relative_tolerance = 1e-8
    hardening_exponent = 2.0
    saturation_resistance = 2.0
    initial_resistance = 0.95
    hardening_modulus = 10.0
    rate_exponent = 0.15
    # reference_strain_rate = 0.05
    activation_energy = 37e3
    gas_constant = 8.314462681 # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = 298

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 0.1
    omega = 13
    alpha = '1 0 0'
    concentration = li_metal_conc
    cref = 0
    block = 'interLayer'
    intBnd = 'blockCeramic_top'
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    enable = true
  [../]

  [./plas2]
    type = ADIsoTropicHyperViscoCreep
    absolute_tolerance = 1e-8
    relative_tolerance = 1e-8

    hardening_exponent = 2.0
    saturation_resistance = 2.0
    initial_resistance = 0.95
    hardening_modulus = 10.0
    rate_exponent = 0.15
    activation_energy = 37e3
    gas_constant = 8.314462681 # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = 298

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 0.1
    block = 'blockMetal'
  [../]
[]

[BCs]
  [./Li_top_y]
    type = ADDirichletBC
    variable = uy
    boundary = 'blockCeramic_bottom'
    value = 0.0
  [../]

  [./Li_top_x]
    type = ADDirichletBC
    boundary = 'blockCeramic_bottom'
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
    type = ADPressure
    boundary = 'blockMetal_top'
    variable = uy
    component = 1
    function = pressure
    use_displaced_mesh = true
  [../]
  # [./li_metal_flux]
  #   type = ScaledCoupledVarNeumannBC
  #   boundary = 'blockMetal_bottom'
  #   variable = li_metal_conc
  #   v = bndliflux
  #   scale = -1.036426e-2
  #   extra_vector_tags = 'ref'
  #   # value = 1.0
  # [../]

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
    block = 'mech_contact_secondary_subdomain'
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

  [./ext_pressure]
    type = SideAverageValue
    variable = stress_yy
    boundary = 'blockMetal_top'
  [../]
  [./matl_ts_min]
    type = MaterialTimeStepPostprocessor
    block = 'interLayer blockMetal'
  [../]
  [./contact_size]
    type = ContactDOFSetSize
    variable = mech_contact_normal_lm
    subdomain = 'mech_contact_secondary_subdomain'
    tolerance = 1e-15
  [../]

  # [./limetal_flux]
  #   type = SideFluxIntegral
  #   variable = li_metal_conc
  #   boundary = 15
  #   diffusivity = diffusivity_Li
  # [../]
[]

[VectorPostprocessors]
  [./metal_flux]
    type = SideValueSampler
    sort_by = x
    boundary = 'blockMetal_bottom'
    variable = 'bndliflux li_metal_conc'
  [../]

  [./current]
    type = SideValueSampler
    sort_by = x
    boundary = 'blockCeramic_top'
    variable = 'bndliflux li_ion_V'
  [../]

  [./interface]
    type = NodalValueSampler
    block = 'mech_contact_secondary_subdomain'
    variable = 'mech_contact_normal_lm thermal_lm'
    sort_by = x
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
  # compute_scaling_once = false
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor'
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1     NONZERO               1e-20               '
  # line_search = basic
  line_search = contact
  l_max_its = 50
  nl_max_its = 35
  nl_abs_tol = 1e-13
  nl_rel_tol = 1e-6
  # start_time = 0.0
  dt = 1.0
  dtmax = 5.0
  dtmin = 1e-5
  end_time = 500.0
  # num_steps = 10
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.5
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 40
    # timestep_limiting_postprocessor = matl_ts_min
  [../]

[] # Executioner


[Outputs]
  [./out]
    type = Exodus
    interval = 2
    file_base = rst/full_model
    sync_times = 1.0
  [../]
  [./csv]
    type = CSV
    file_base = csv/full_model
  [../]
  [./checkpoint]
    type = Checkpoint
    file_base = check/full_model
    num_files = 4
    interval = 10
    sync_times = 1.0
  [../]

[] # Outputs

[Debug]
  show_var_residual_norms = true
  show_material_props = true
[]
