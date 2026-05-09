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
    new_block_name = 'LiLayer AgC'
    input = 'outflow_area'
  []
  [interface]
    type = SideSetsBetweenSubdomainsGenerator
    new_boundary = 'interface'
    primary_block = 'AgC'
    paired_block = 'LiLayer'
    input = rename
  []
[]

[Problem]
  coord_type = RZ
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

[GlobalParams]
  PorousFlowDictator = dictator
  biot_coefficient = 1.0
  displacements = 'disp_r disp_z'
  gravity = '0 0 0'
[]

[Variables]
  [porepressure]
  []
  [disp_r]
  []
  [disp_z]
  []
[]

[Kernels]
  [time_derivative]
    type = PorousFlowMassTimeDerivative
    variable = porepressure
  []
  [flux]
    type = PorousFlowAdvectiveFlux
    variable = porepressure
    gravity = '0 0 0'
    use_displaced_mesh = false
  []
  [vol_strain_rate]
    type = PorousFlowMassVolumetricExpansion
    fluid_component = 0
    variable = porepressure
  []
  [grad_stress_r]
    type = StressDivergenceRZTensors
    # temperature = T
    variable = disp_r
    # eigenstrain_names = thermal_contribution
    use_displaced_mesh = false
    component = 0
  []

  [effective_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
  []

  [poro_r]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_r
    use_displaced_mesh = false
    component = 0
  []
  [grad_stress_z]
    type = StressDivergenceRZTensors
    # temperature = T
    variable = disp_z
    # eigenstrain_names = thermal_contribution
    use_displaced_mesh = false
    component = 1
  []
  [poro_z]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_z
    use_displaced_mesh = false
    component = 1
  []
[]

[BCs]
  [constant_injection_porepressure]
    type = DirichletBC
    variable = porepressure
    value = 1e5
    boundary = top
  []
  [constant_injection_flux]
    type = PorousFlowSink
    variable = porepressure
    flux_function = "if (t <= 500,if (t <10, -66.41344e-9*3*t,-66.41344e-9*30), 66.41344e-9*30)"
    # flux_function = -664.13e-9
    boundary = bottom
    fluid_phase = 0
    use_relperm = true
  []
  [pinned_left_r]
    type = DirichletBC
    variable = disp_r
    value = 0
    boundary = 'left'
  []
  [pinned_bottom_z]
    type = DirichletBC
    variable = disp_z
    value = 0
    boundary = 'bottom'
  []
[]

[AuxVariables]
  [sat]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_rr]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_tt]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_zz]
    family = MONOMIAL
    order = CONSTANT
  []
  [porosity]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [saturation]
    type = PorousFlowPropertyAux
    variable = sat
    property = saturation
  []
  [stress_rr_aux]
    type = RankTwoAux
    variable = stress_rr
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
  []
  [stress_tt]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_tt
    index_i = 2
    index_j = 2
  []
  [stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 1
    index_j = 1
  []
  [porosity]
    type = PorousFlowPropertyAux
    variable = porosity
    property = porosity
    execute_on = timestep_end
  []
[]

[Modules]
  [FluidProperties]
    [the_simple_fluid]
      type = SimpleFluidProperties
      bulk_modulus = 11E9
      viscosity = 1.0e2
      density0 = 534
      thermal_expansion = 0.0000118
      cp = 4194
      cv = 4186
      porepressure_coefficient = 0
    []
  []
[]

[Materials]
  [porosity_Agc]
    type = PorousFlowPorosityConst
    porosity = 0.5
    chemical = false
    block = AgC
    use_displaced_mesh = false
    # initial_mineral_concentrations = initial_and_reference_conc
    # reference_chemistry = initial_and_reference_conc
  []

  [porosity_LiLayer]
    type = PorousFlowPorosityConst
    porosity = 0.95
    chemical = false
    block = LiLayer
    # initial_mineral_concentrations = initial_and_reference_conc
    # reference_chemistry = initial_and_reference_conc
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
  [saturation_calculator]
    type = PorousFlow1PhaseP
    porepressure = porepressure
    capillary_pressure = pc
  []

  [massfrac]
    type = PorousFlowMassFraction
  []
  [simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = the_simple_fluid
    phase = 0
  []
  [relperm]
    type = PorousFlowRelativePermeabilityCorey
    n = 3
    s_res = 0.1
    sum_s_res = 0.1
    phase = 0
  []
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 5E9
    poissons_ratio = 0.25
  []
  [strain]
    type = ComputeAxisymmetricRZSmallStrain
  []
  [stress]
    type = ComputeLinearElasticStress
  []
  [volumetric_strain]
    type = PorousFlowVolumetricStrain
  []
  [effective_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
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
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_force_iteration'
    petsc_options_value = ' lu       mumps 1'
  []
[]

[Executioner]
  type = Transient
  automatic_scaling = true
  # compute_scaling_once = false
  solve_type = Newton
  # end_time = 36000
  dt = 10
  nl_max_its = 100
  nl_abs_tol = 1E-8
  nl_rel_tol = 1E-4
  resid_vs_jac_scaling_param = 0.5
  dtmax = 10
  end_time = 1000
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0
    optimal_iterations = 100
    cutback_factor = 0.5
    growth_factor = 1.1
  []

[]

[Outputs]
  exodus = true
  csv = true
  sync_times = '500'
  sync_only = false
  file_base = xxx2
[]
