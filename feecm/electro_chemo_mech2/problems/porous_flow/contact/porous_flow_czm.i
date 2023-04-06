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
[Mesh]
  [file]
    type = FileMeshGenerator
    file = test_czm.e
  []
  [break]
    type = BreakMeshByBlockGenerator
    input = file
    split_interface = true
  []
[]

[Problem]
  # coord_type = RZ
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  acceptable_iterations = 2
  group_variables = 'disp_x disp_y porepressure; li_ion_V'
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  include_gap = false
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [li_ion_V]
    block = 'Anode  Electrolyte'
  []
[]

# [Constraints]
  # [thermal_constraint]
  #   type = GapDisplacementConductanceConstraint
  #   variable = thermal_lm
  #   secondary_variable = li_ion_V
  #   primary_boundary = 'e_top'
  #   secondary_subdomain = 'secondary'
  #   secondary_boundary = 'p_bottom'
  #   primary_subdomain = 'primary'
  #   # k_function = gapk
  #   k = 1e2
  #   use_displaced_mesh = true
  #   compute_lm_residuals = true
  #   extra_vector_tags = 'ref'
  #   include_gap = true
  #   include_equilibrium_potential = false
  #   R = 8.3145
  #   faraday = 96485.3329
  #   temperature = 298
  #   surfaceType = SECONDARY
  #   include_concentration = false
  #
  # []
  # [normal_lm_x]
  #   type = EqualValueConstraint
  #   variable = normal_lm_x
  #   primary_variable = disp_r
  #   secondary_variable = disp_r
  #   primary_boundary = 'e_top'
  #   secondary_subdomain = 'secondary'
  #   secondary_boundary = 'p_bottom'
  #   primary_subdomain = 'primary'
  #   use_displaced_mesh = false
  #   # delta = 1.0
  #   # delta = 0.1
  # []
  #
  # [normal_lm_y]
  #   type = EqualValueConstraint
  #   variable = normal_lm_y
  #   primary_variable = disp_z
  #   secondary_variable = disp_z
  #   primary_boundary = 'e_top'
  #   secondary_subdomain = 'secondary'
  #   secondary_boundary = 'p_bottom'
  #   primary_subdomain = 'primary'
  #   use_displaced_mesh = false
  #   # delta = 0.5
  # []
  # [pp_lm]
  #   type = EqualValueConstraint
  #   variable = pp_lm
  #   primary_variable = porepressure
  #   secondary_variable = porepressure
  #   # component = 1
  #   primary_boundary = 'e_top'
  #   secondary_subdomain = 'secondary'
  #   secondary_boundary = 'p_bottom'
  #   primary_subdomain = 'primary'
  #   use_displaced_mesh = false
  # []
#   [normal_lm]
#     type = NormalNodalLMMechanicalContact
#     secondary = 'p_bottom'
#     primary = 'e_top'
#     variable = normal_lm
#     primary_variable = disp_r
#     normal_smoothing_distance = 0.2
#     disp_y = disp_z
#     ncp_function_type = min
#     tangential_tolerance = 0.2
#     # c = 1e-3
#   []
#   [normal_r]
#     type = NormalMortarMechanicalContact
#     primary_boundary = 'e_top'
#     secondary_boundary = 'p_bottom'
#     primary_subdomain = 'primary'
#     secondary_subdomain = 'secondary'
#     variable = normal_lm
#     secondary_variable = disp_r
#     component = x
#     use_displaced_mesh = true
#     compute_lm_residuals = false
#     # extra_vector_tags = 'ref'
#   []
#   [normal_z]
#     type = NormalMortarMechanicalContact
#     secondary_boundary = 'e_top'
#     primary_boundary = 'p_bottom'
#     primary_subdomain = 'primary'
#     secondary_subdomain = 'secondary'
#     variable = normal_lm
#     secondary_variable = disp_z
#     component = y
#     use_displaced_mesh = true
#     compute_lm_residuals = false
#     # extra_vector_tags = 'ref'
#   []
# []
[ICs]
  [porepressure]
    type = ConstantIC
    variable = porepressure
    block = 'Anode '
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
#     # block = 'Anode  Electrolyte'
#     # block = 'Anode  Electrolyte'
#     extra_vector_tags = 'ref'
#     # use_finite_deform_jacobian = true
#   []
# []

[Kernels]
  # [TensorMechanics]
  # []
  [li_ion_V]
    type = ADHeatConduction
    block = 'Anode  Electrolyte'
    variable = li_ion_V
    use_displaced_mesh = false
  []
  [time_derivative]
    type = PorousFlowMassTimeDerivative
    variable = porepressure
    fluid_component = 0
    block = 'Anode '
  []
  [flux]
    type = PorousFlowAdvectiveFlux
    variable = porepressure
    gravity = '0 0 0'
    fluid_component = 0
    block = 'Anode '
  []
  [vol_strain_rate_water]
    type = PorousFlowMassVolumetricExpansion
    fluid_component = 0
    variable = porepressure
    block = 'Anode '
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
    block = 'Anode '
    # block = 'Anode '
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
    block = 'Anode '
    # block = 'Anode '
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
  [top_mech_pressure]
    type = DirichletBC
    variable = porepressure
    value = '1e5'
    boundary = p_top
  []
  [pinned_top_bottom_z]
    type = DirichletBC
    variable = disp_z
    value = 0
    boundary = 'e_bottom'
  []
  # [constant_injection_porepressure]
  #   type = DirichletBC
  #   variable = porepressure
  #   value = 1e5
  #   boundary = p_top
  # []
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
[]

[Materials]
  [thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1.0e-2
    # block = 'blockCeramic'
  []

  [porosity_Anode]
    type = PorousFlowPorosity
    ensure_positive = true
    mechanical = true
    porosity_zero = 0.5
    chemical = false
    fluid = true
    block = Anode
    solid_bulk = 100e9
    # initial_mineral_concentrations = initial_and_reference_conc
    # reference_chemistry = initial_and_reference_conc
  []

  [porosity_]
    type = PorousFlowPorosity
    porosity_zero = 0.95
    mechanical = true
    ensure_positive = true
    fluid = true
    chemical = false
    block =
    solid_bulk = 5e9
    # initial_mineral_concentrations = initial_and_reference_conc
    # reference_chemistry = initial_and_reference_conc
  []
  # [porosity_Electrolyte]
  #   type = PorousFlowPorosityConst
  #   block = 'Electrolyte secondary primary'
  #   porosity = 1.0
  # []

  [permeability_Anode]
    type = PorousFlowPermeabilityKozenyCarman
    block = Anode
    k_anisotropy = '1E-25 0 0   0 1e-1 0   0 0 1E-25'
    f = 0.01
    d = 50e-9
    poroperm_function = kozeny_carman_fd2
    m = 3
    n = 2
  []
  [permeability_]
    type = PorousFlowPermeabilityKozenyCarman
    block =
    k_anisotropy = '1 0 0   0 1 0   0 0 1'
    f = 0.01
    d = 50e-9
    poroperm_function = kozeny_carman_fd2
    m = 3
    n = 2
  []

  # [permeability_Electrolyte]
  #   type = PorousFlowPermeabilityConst
  #   block = 'Electrolyte secondary primary'
  #   permeability = '0 0 0   0 0 0   0 0 0'
  # []
  [saturation_calculator]
    type = PorousFlow1PhaseP
    porepressure = porepressure
    capillary_pressure = pc
    block = 'Anode '
  []
  [temperature]
    type = PorousFlowTemperature
    temperature = 293
    block = 'Anode '
  []
  [massfrac]
    type = PorousFlowMassFraction
    block = 'Anode '
  []
  [simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = the_simple_fluid
    phase = 0
    block = 'Anode '
  []
  [relperm]
    type = PorousFlowRelativePermeabilityCorey
    n = 3
    s_res = 0.001
    sum_s_res = 0.001
    phase = 0
    block = 'Anode '
  []
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 100E9
    poissons_ratio = 0.25
    # block = 'Anode  Electrolyte'
  []
  [strain]
    type = ComputeSmallStrain
    # eigenstrain_names = 'thermal_contribution initial_stress'
  []
  [stress]
    type = ComputeLinearElasticStress
    # block = 'Anode  Electrolyte'
  []
  [effective_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
    block = 'Anode '
  []
  [volumetric_strain]
    type = PorousFlowVolumetricStrain
    block = 'Anode '
  []
  [biot_modulus]
    type = PorousFlowConstantBiotModulus
    solid_bulk_compliance = 2E-9
    fluid_bulk_modulus = 11E9
    block = 'Anode '
  []
  [density_Anode_Li]
    type = GenericConstantMaterial
    block = 'Anode '
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
  nl_rel_tol = 1E-5
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
  file_base = xxx3_new
[]

[VectorPostprocessors]
  [porous_bottom]
    type = SideValueSampler
    sort_by = x
    boundary = p_bottom
    variable = 'bndliflux li_ion_V li_ion_flux_y nodal_outflow'
  []
  [Electrolyte_top]
    type = SideValueSampler
    sort_by = x
    boundary = p_bottom
    variable = 'bndliflux li_ion_V li_ion_flux_y'
  []
[]

[Postprocessors]
  [bottom_current]
    type = ADSideDiffusiveFluxIntegral
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
  [contact]
    type = ContactDOFSetSize
    variable = normal_lm
    subdomain = secondary
  []
[]
