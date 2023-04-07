# Darcy flow with heat advection and conduction, and elasticity
# Dissolution of Ag in Li is the problem we are trying to solve here
# Model consists of a single block of porous material
#   1) Single phase fluid (Li/ LiAg)
#   2) Ag dissolves in the fluid
# If we want to model Ag precipitation in the fluid then perhaps we need 2 phases
# Li metal viscosity is set to 1e4
# Conversion from current density to mass flux = i A/m^2 / F(C/mol) * rho (kg/m^2) * omega (m^3/mol) = kg/m^2s
#                                               = current_density / Faradays' * molar volume * density
#                                               = 1 A/m^2 / 96485.2239 (C/mol) * 1.2e-5 (m^3/mol) * 534 (kg/m^3)
#                                               = 66.143144e-9 kg/m^2s
# Permeability in the Porous layer is Carman Kozeny formulation of type 1 where k_ij = A k_ij^0 * phi^n/ (1-phi)^m
#                                                                   A = f d^2 -> d is a grain size or particle size
# [Mesh]
#   [file]
#     type = FileMeshGenerator
#     file = test2.e
#   []
#   [secondary_boundary_block]
#     type = LowerDBlockFromSidesetGenerator
#     input = file
#     sidesets = 'p_bottom'
#     new_block_name = 'secondary'
#   []
#   [primary_boundary_block]
#     type = LowerDBlockFromSidesetGenerator
#     input = secondary_boundary_block
#     sidesets = 'e_top'
#     new_block_name = 'primary'
#   []
# []
[Mesh]
  [mesh]
    # type = FileMeshGenerator
    # file = test_mesh.e
    type = GeneratedMeshGenerator
    dim = 2
    nx = 60
    xmin = 0.0
    xmax = 100e-6
    bias_x = 1
    ny = 30
    ymin = 0
    ymax = 30e-6
    boundary_name_prefix = 'p'
  []

  [AgC]
    type = SubdomainBoundingBoxGenerator
    block_id = 1
    bottom_left = '0.0 0.0 0'
    top_right = '100e-6 15.01e-6 0'
    input = mesh
  []
  [injection_area]
    type = ParsedGenerateSideset
    # combinatorial_geometry = 'x*x+y*y<1.01'
    combinatorial_geometry = 'y<0.0001e-6'
    included_subdomain_ids = 1
    new_sideset_name = 'injection_area'
    input = 'AgC'
  []
  [outflow_area]
    type = ParsedGenerateSideset
    combinatorial_geometry = 'y>29.99e-6'
    included_subdomain_ids = 1
    new_sideset_name = 'outflow_area'
    input = injection_area
  []
  [rename]
    type = RenameBlockGenerator
    old_block_id = '0 1'
    new_block_id = '100 101'
    input = 'outflow_area'
  []
  [rename2]
    type = RenameBlockGenerator
    old_block_id = '100 101'
    input = rename
    new_block_name = 'LiLayer AgC'
  []

  [electrolyte]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 45
    ny = 15
    xmin = 0.0
    xmax = 100e-6
    ymin = -20e-6
    ymax = 0.0
    boundary_name_prefix = 'e'
    boundary_id_offset = 100
  []
  [rename_electrolyte]
    type = RenameBlockGenerator
    old_block_id = 0
    new_block_id = 200
    input = electrolyte
  []
  [full_mesh]
    type = MeshCollectionGenerator
    inputs = 'rename2 rename_electrolyte'
  []
  # [full_mesh]
  #   type = StitchedMeshGenerator
  #   inputs = 'rename2 rename_electrolyte'
  #   stitch_boundaries_pairs = 'p_bottom e_top'
  #   clear_stitched_boundary_ids = false
  # []
  [rename_electrolyte2]
    type = RenameBlockGenerator
    old_block_id = 200
    new_block_name = 'electrolyte'
    input = full_mesh
  []
  [secondary_boundary_block]
    type = LowerDBlockFromSidesetGenerator
    input = rename_electrolyte2
    sidesets = 'p_bottom'
    new_block_name = 'secondary'
  []
  [primary_boundary_block]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_boundary_block
    sidesets = 'e_top'
    new_block_name = 'primary'
  []
  # [interface]
  #   type = BreakMeshByBlockGenerator
  #   input = rename_electrolyte2
  #   block_pairs = '101 200'
  #   add_interface_on_two_sides = true
  #   split_interface = true
  #   show_info = true
  # []
  # [common_boundary]
  #   type = RenameBoundaryGenerator
  #   input = rename_electrolyte2
  #   old_boundary = 'p_bottom e_top'
  #   new_boundary = 'interface interface'
  #   show_info = true
  # []

[]
# [ICs]
#   [porepressure]
#     type = ConstantIC
#     value = -1e-6
#     variable = porepressure
#   []
# []
[Problem]
  # coord_type = RZ
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  acceptable_iterations = 2
  group_variables = 'disp_r disp_z normal_lm_x normal_lm_y porepressure; li_ion_V thermal_lm'

[]

[GlobalParams]
  PorousFlowDictator = dictator
  biot_coefficient = 1.0
  displacements = 'disp_r disp_z'
  gravity = '0 0 0'
[]

[AuxVariables]
  [nodal_outflow]
  []
  [sat]
    family = MONOMIAL
    order = CONSTANT
    block = 'AgC LiLayer'
  []
  [permeability]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_xx]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_yy]
    family = MONOMIAL
    order = CONSTANT
  []
  [hydrostatic_stress]
    family = MONOMIAL
    order = CONSTANT
  []
  [vonmises_stress]
    family = MONOMIAL
    order = CONSTANT
  []
  [bndliflux]
    order = CONSTANT
    family = MONOMIAL
    block = 'AgC LiLayer electrolyte'
    # block = 'blockCeramic interLayer blockMetal blockQuartz'
  []
  [li_ion_flux_x]
    order = FIRST
    family = MONOMIAL
    block = 'AgC LiLayer electrolyte'
    # block = 'blockCeramic interLayer blockMetal blockQuartz'
  []
  [li_ion_flux_y]
    order = FIRST
    family = MONOMIAL
    block = 'AgC LiLayer electrolyte'
    # block = 'blockCeramic interLayer blockMetal blockQuartz'
  []
  [li_ion_flux_z]
    order = FIRST
    family = MONOMIAL
    block = 'AgC LiLayer electrolyte'
    # block = 'blockCeramic interLayer blockMetal blockQuartz'
  []
  [darcy_vel_x]
    type = MooseVariableConstMonomial
  []
  [darcy_vel_y]
    type = MooseVariableConstMonomial
  []
  [darcy_vel_z]
    type = MooseVariableConstMonomial
  []
  [effective_fluid_pressure]
    family = MONOMIAL
    order = CONSTANT
    block = 'AgC LiLayer'
  []

[]

[AuxKernels]
  [effective_fluid_pressure]
    type = ParsedAux
    args = 'porepressure sat'
    function = 'porepressure * sat'
    variable = effective_fluid_pressure
    block = 'AgC LiLayer'
  []
  [bnd_li_flux]
    type = DiffusionFluxNormalToBoundaryAux
    variable = bndliflux
    boundary = 'e_top p_bottom'
    diffusion_variable = li_ion_V
    diffusivity = thermal_conductivity
    block = 'AgC LiLayer electrolyte'
  []
  [li_ion_flux_x]
    type = ADDiffusionFluxAux
    variable = li_ion_flux_x
    component = x
    diffusion_variable = li_ion_V
    diffusivity = thermal_conductivity
    block = 'AgC LiLayer electrolyte'
    # block = 'blockCeramic blockMetal interLayer blockQuartz'
  []

  [li_ion_flux_y]
    type = ADDiffusionFluxAux
    variable = li_ion_flux_y
    component = y
    diffusion_variable = li_ion_V
    diffusivity = thermal_conductivity
    block = 'AgC LiLayer electrolyte'
    # block = 'blockCeramic blockMetal interLayer blockQuartz'
  []
  [li_ion_flux_z]
    type = ADDiffusionFluxAux
    variable = li_ion_flux_z
    component = z
    diffusion_variable = li_ion_V
    diffusivity = thermal_conductivity
    block = 'AgC LiLayer electrolyte'
    # block = 'blockCeramic blockMetal interLayer blockQuartz'
  []

  [saturation]
    type = PorousFlowPropertyAux
    variable = sat
    property = saturation
    block = 'AgC LiLayer'
  []
  [permeability]
    type = PorousFlowPropertyAux
    property = permeability
    column = 1
    row = 1
    variable = permeability
  []
  [stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  []
  [stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  []
  [hydrostatic_stress]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    scalar_type = Hydrostatic
    variable = hydrostatic_stress
  []
  [vonmises_stress]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    scalar_type = VonMisesStress
    variable = vonmises_stress
  []
  [PorousFlowActionBase_Darcy_x_Aux]
    type = PorousFlowDarcyVelocityComponent
    PorousFlowDictator = dictator
    component = x
    execute_on = TIMESTEP_END
    variable = darcy_vel_x
  []
  [PorousFlowActionBase_Darcy_y_Aux]
    type = PorousFlowDarcyVelocityComponent
    PorousFlowDictator = dictator
    component = y
    execute_on = TIMESTEP_END
    variable = darcy_vel_y
  []
  [PorousFlowActionBase_Darcy_z_Aux]
    type = PorousFlowDarcyVelocityComponent
    PorousFlowDictator = dictator
    component = z
    execute_on = TIMESTEP_END
    variable = darcy_vel_z
  []

[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'porepressure disp_r disp_z'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
  [pc]
    type = PorousFlowCapillaryPressureVG
    alpha = 1E-6
    m = 0.6
  []
[]

[Modules]
  [FluidProperties]
    [the_simple_fluid]
      type = SimpleFluidProperties
      bulk_modulus = 11E9
      viscosity = 1.0e2
      density0 = 534
      thermal_expansion = 0
      cp = 4194
      cv = 4186
      porepressure_coefficient = 0
    []
  []
[]

[Variables]
  [disp_r]
  []
  [disp_z]
  []
  [normal_lm_x]
    block = 'secondary'
  []
  [normal_lm_y]
    block = 'secondary'
  []

  # [normal_lm]
  #   block = 'secondary'
  # []
  [thermal_lm]
    block = 'secondary'
  []
  [li_ion_V]
    block = 'AgC LiLayer electrolyte'
  []
  [porepressure]
  []
[]

[Constraints]
  [thermal_constraint]
    type = GapDisplacementConductanceConstraint
    variable = thermal_lm
    secondary_variable = li_ion_V
    primary_boundary = 'e_top'
    secondary_subdomain = 'secondary'
    secondary_boundary = 'p_bottom'
    primary_subdomain = 'primary'
    # k_function = gapk
    k = 1e2
    use_displaced_mesh = false
    compute_lm_residuals = true
    extra_vector_tags = 'ref'
    include_gap = false
    include_equilibrium_potential = false
    R = 8.3145
    faraday = 96485.3329
    temperature = 298
    surfaceType = SECONDARY
    include_concentration = false

  []
  [normal_lm_x]
    type = EqualValueConstraint
    variable = normal_lm_x
    primary_variable = disp_r
    secondary_variable = disp_r
    primary_boundary = 'e_top'
    secondary_subdomain = 'secondary'
    secondary_boundary = 'p_bottom'
    primary_subdomain = 'primary'
    use_displaced_mesh = true
    # delta = 1.0
    # delta = 0.1
  []

  [normal_lm_y]
    type = EqualValueConstraint
    variable = normal_lm_y
    primary_variable = disp_z
    secondary_variable = disp_z
    primary_boundary = 'e_top'
    secondary_subdomain = 'secondary'
    secondary_boundary = 'p_bottom'
    primary_subdomain = 'primary'
    use_displaced_mesh = true
    # delta = 0.5
  []
  # [normal_lm]
  #   type = NormalNodalLMMechanicalContact
  #   secondary = 'p_bottom'
  #   primary = 'e_top'
  #   variable = normal_lm
  #   primary_variable = disp_r
  #   normal_smoothing_distance = 0.2
  #   disp_y = disp_z
  #   ncp_function_type = min
  #   tangential_tolerance = 0.2
  #   # c = 1e-3
  # []
  # [normal_r]
  #   type = NormalMortarMechanicalContact
  #   primary_boundary = 'e_top'
  #   secondary_boundary = 'p_bottom'
  #   primary_subdomain = 'primary'
  #   secondary_subdomain = 'secondary'
  #   variable = normal_lm
  #   secondary_variable = disp_r
  #   component = x
  #   use_displaced_mesh = true
  #   compute_lm_residuals = false
  #   # extra_vector_tags = 'ref'
  # []
  # [normal_z]
  #   type = NormalMortarMechanicalContact
  #   secondary_boundary = 'e_top'
  #   primary_boundary = 'p_bottom'
  #   primary_subdomain = 'primary'
  #   secondary_subdomain = 'secondary'
  #   variable = normal_lm
  #   secondary_variable = disp_z
  #   component = y
  #   use_displaced_mesh = true
  #   compute_lm_residuals = false
  #   # extra_vector_tags = 'ref'
  # []
[]
[ICs]
  [porepressure]
    type = ConstantIC
    variable = porepressure
    block = 'AgC LiLayer'
    value = 1e5
  []
[]
# [Modules/TensorMechanics/Master]
#   [all]
#     add_variables = true
#     strain = SMALL
#     volumetric_locking_correction = false
#     generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress hydrostatic_stress'
#     use_automatic_differentiation = false
#     # block = 'AgC LiLayer electrolyte'
#     # block = 'AgC LiLayer electrolyte'
#     extra_vector_tags = 'ref'
#     # use_finite_deform_jacobian = true
#   []
# []

[Kernels]
  # [TensorMechanics]
  # []
  [li_ion_V]
    type = ADHeatConduction
    block = 'AgC LiLayer electrolyte'
    variable = li_ion_V
    use_displaced_mesh = false
  []
  [time_derivative]
    type = PorousFlowMassTimeDerivative
    variable = porepressure
    fluid_component = 0
  []
  [flux]
    type = PorousFlowAdvectiveFlux
    variable = porepressure
    gravity = '0 0 0'
    fluid_component = 0
  []
  [vol_strain_rate_water]
    type = PorousFlowMassVolumetricExpansion
    fluid_component = 0
    variable = porepressure
    block = 'AgC LiLayer'
  []
  [grad_stress_r]
    type = StressDivergenceTensors
    # temperature = T
    variable = disp_r
    # eigenstrain_names = thermal_contribution
    use_displaced_mesh = false
    component = 0
    extra_vector_tags = 'ref'
  []
  [poro_r]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_r
    use_displaced_mesh = false
    component = 0
    # block = 'AgC LiLayer'
  []
  [grad_stress_z]
    type = StressDivergenceTensors
    # temperature = T
    variable = disp_z
    # eigenstrain_names = thermal_contribution
    use_displaced_mesh = false
    component = 1
    extra_vector_tags = 'ref'
  []
  [poro_z]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_z
    use_displaced_mesh = false
    component = 1
    # block = 'AgC LiLayer'
  []
[]

[BCs]
  [constant_current_density]
    type = ADFunctionNeumannBC
    boundary = 'e_bottom'
    variable = li_ion_V
    function = 'if (t <= 10, 1.0*t, 10.0 )'
    extra_vector_tags = 'ref'
  []
  [constant_Voltage]
    type = DirichletBC
    boundary = 'p_top'
    variable = li_ion_V
    value = 0.0
  []
  [pinned_top_bottom_r]
    type = DirichletBC
    variable = disp_r
    value = 0
    boundary = 'p_left e_left p_right e_right'
  []
  [pinned_top_bottom_z]
    type = DirichletBC
    variable = disp_z
    value = 0
    boundary = 'p_top e_bottom'
  []
  [constant_injection_porepressure]
    type = DirichletBC
    variable = porepressure
    value = 1e5
    boundary = p_top
  []
  [constant_injection_flux_scaled]
    type = PorousFlowSinkScaledCoupledVar
    v = bndliflux
    variable = porepressure
    scale = 66.41344e-9
    flux_function = '1.0'
    boundary = p_bottom
    fluid_phase = 0
    # use_relperm = true
    extra_vector_tags = 'ref'
    save_in = nodal_outflow
  []
  [constant_injection_flux_scaled_electrolyte]
    type = PorousFlowSink
    # v = bndliflux
    variable = porepressure
    # scale = 66.41344e-9
    flux_function = '0.0'
    boundary = e_top
    fluid_phase = 0
    # use_relperm = true
    extra_vector_tags = 'ref'
    save_in = nodal_outflow
  []
[]

[Materials]
  [thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1.0e-2
    # block = 'blockCeramic'
  []

  [porosity_Agc]
    type = PorousFlowPorosity
    ensure_positive = true
    mechanical = true
    porosity_zero = 0.5
    chemical = false
    fluid = true
    block = AgC
    solid_bulk = 100e9
    # initial_mineral_concentrations = initial_and_reference_conc
    # reference_chemistry = initial_and_reference_conc
  []

  [porosity_LiLayer]
    type = PorousFlowPorosity
    porosity_zero = 0.95
    mechanical = true
    ensure_positive = true
    fluid = true
    chemical = false
    block = LiLayer
    solid_bulk = 5e9
    # initial_mineral_concentrations = initial_and_reference_conc
    # reference_chemistry = initial_and_reference_conc
  []
  [porosity_electrolyte]
    type = PorousFlowPorosityConst
    block = 'electrolyte secondary primary'
    porosity = 1.0
  []

  [permeability_AgC]
    type = PorousFlowPermeabilityKozenyCarman
    block = AgC
    k_anisotropy = '1E-25 0 0   0 1 0   0 0 1E-25'
    f = 0.01
    d = 50e-9
    poroperm_function = kozeny_carman_fd2
    m = 3
    n = 2
  []
  [permeability_LiLayer]
    type = PorousFlowPermeabilityKozenyCarman
    block = LiLayer
    k_anisotropy = '1 0 0   0 1 0   0 0 1'
    f = 0.01
    d = 50e-9
    poroperm_function = kozeny_carman_fd2
    m = 3
    n = 2
  []

  [permeability_electrolyte]
    type = PorousFlowPermeabilityConst
    block = 'electrolyte secondary primary'
    permeability = '0 0 0   0 0 0   0 0 0'
  []
  [saturation_calculator]
    type = PorousFlow1PhaseP
    porepressure = porepressure
    capillary_pressure = pc
    # block = 'AgC LiLayer'
  []
  [temperature]
    type = PorousFlowTemperature
    temperature = 293
    # block = 'AgC LiLayer'
  []
  [massfrac]
    type = PorousFlowMassFraction
    # block = 'AgC LiLayer'
  []
  [simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = the_simple_fluid
    phase = 0
    # block = 'AgC LiLayer'
  []
  [relperm]
    type = PorousFlowRelativePermeabilityCorey
    n = 3
    s_res = 0.001
    sum_s_res = 0.001
    phase = 0
    # block = 'AgC LiLayer'
  []
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 100E9
    poissons_ratio = 0.25
    # block = 'AgC LiLayer electrolyte'
  []
  [strain]
    type = ComputeSmallStrain
    # eigenstrain_names = 'thermal_contribution initial_stress'
  []
  [stress]
    type = ComputeLinearElasticStress
    # block = 'AgC LiLayer electrolyte'
  []
  [effective_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
    # block = 'AgC LiLayer'
  []
  [volumetric_strain]
    type = PorousFlowVolumetricStrain
    # block = 'AgC LiLayer'
  []
  [biot_modulus]
    type = PorousFlowConstantBiotModulus
    solid_bulk_compliance = 2E-9
    fluid_bulk_modulus = 11E9
    # block = 'AgC LiLayer'
  []
  [density_AgC_Li]
    type = GenericConstantMaterial
    # block = 'AgC LiLayer'
    prop_names = density
    prop_values = 534
  []

[]

[Preconditioning]
  active = preferred_but_might_not_be_installed
  [basic]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = ' asm      lu           NONZERO                   2'
  []
  [preferred_but_might_not_be_installed]
    type = SMP
    full = true
    # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_force_iteration'
    # petsc_options_value = ' lu       mumps 1'
  []
[]

[Executioner]
  type = Transient
  automatic_scaling = true
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor'
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration '
                        '-pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1     NONZERO               1e-20               '

  # compute_scaling_once = false
  solve_type = Newton
  # end_time = 36000
  dt = 10
  nl_max_its = 100
  nl_abs_tol = 1E-12
  nl_rel_tol = 1E-6
  resid_vs_jac_scaling_param = 0.5
  dtmax = 10
  end_time = 500
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0
    optimal_iterations = 100
    cutback_factor = 0.5
    growth_factor = 1.1
  []
  verbose = true

[]

[Outputs]
  exodus = true
  csv = true
  execute_on = 'INITIAL TIMESTEP_END'
  # sync_times = '500'
  # sync_only = false
  file_base = xxx3
[]

[VectorPostprocessors]
  [porous_bottom]
    type = SideValueSampler
    sort_by = x
    boundary = p_bottom
    variable = 'bndliflux li_ion_V li_ion_flux_y nodal_outflow'
  []
  [electrolyte_top]
    type = SideValueSampler
    sort_by = x
    boundary = p_bottom
    variable = 'bndliflux li_ion_V li_ion_flux_y'
  []
[]

[Postprocessors]
  [bottom_current]
    type = ADSideFluxIntegral
    variable = li_ion_V
    boundary = 'e_bottom'
    diffusivity = thermal_conductivity
  []
  [stress_zz_bottom]
    type = SideAverageValue
    variable = stress_yy
    boundary = p_bottom
  []
  [saturation]
    type = SideAverageValue
    variable = sat
    boundary = p_bottom
  []
  [porepressure_bottom]
    type = SideAverageValue
    variable = porepressure
    boundary = p_bottom
  []
[]
