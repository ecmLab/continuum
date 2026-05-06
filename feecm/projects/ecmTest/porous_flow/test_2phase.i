# Two phase, temperature-independent, with mechanics, radial with fine mesh, constant injection of Li into a porous interlayer
# species=0 is Li
# species=1 is gas
# phase=0 is liquid, and since massfrac_ph0_sp0 = 1, this is all Li
# phase=1 is gas, and since massfrac_ph1_sp0 = 0, this is all gas
#
# The mesh used below has very high resolution, so the simulation takes a long time to complete.
# Some suggested meshes of different resolution:
# nx=50, bias_x=1.2
# nx=100, bias_x=1.1
# nx=200, bias_x=1.05
# nx=400, bias_x=1.02
# nx=1000, bias_x=1.01
# nx=2000, bias_x=1.003
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 500
  bias_x = 1.003
  xmin = 0.1e-6
  xmax = 10e-6
  ny = 1
  ymin = 0
  ymax = 0.1e-6
[]

[Problem]
  coord_type = RZ
[]

[GlobalParams]
  displacements = 'disp_r disp_z'
  PorousFlowDictator = dictator
  gravity = '0 0 0'
  biot_coefficient = 1.0
[]

[Variables]
  [pli]
    initial_condition = 0.0
    # initial_condition = 18.3e6
  []
  [gas_saturation]
    initial_condition = 0.0
  []
  [temp]
    initial_condition = 293
  []
  [disp_r]
  []
  # [disp_z]
  # []

[]

[AuxVariables]
  [rate]
  []
  [disp_z]
  []
  [massfrac_ph0_sp0]
    initial_condition = 1 # all Li in phase=0
  []
  [massfrac_ph1_sp0]
    initial_condition = 0 # no Li in phase=1
  []
  [pgas]
    family = MONOMIAL
    order = FIRST
  []
  [sLi]
    family = MONOMIAL
    order = FIRST
  []
  [stress_rr]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_tt]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_zz]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [mass_Li_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pli
  []
  [flux_Li]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    use_displaced_mesh = false
    variable = pli
  []
  [mass_gas_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = gas_saturation
  []
  [flux_gas]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    use_displaced_mesh = false
    variable = gas_saturation
  []
  [energy_dot]
    type = PorousFlowEnergyTimeDerivative
    variable = temp
  []
  [advection]
    type = PorousFlowHeatAdvection
    use_displaced_mesh = false
    variable = temp
  []
  [conduction]
    type = PorousFlowExponentialDecay
    use_displaced_mesh = false
    variable = temp
    reference = 293
    rate = rate
  []
  [grad_stress_r]
    type = StressDivergenceRZTensors
    temperature = temp
    eigenstrain_names = thermal_contribution
    variable = disp_r
    use_displaced_mesh = false
    component = 0
  []
  # [grad_stress_z]
  #   type = StressDivergenceRZTensors
  #   temperature = temp
  #   eigenstrain_names = thermal_contribution
  #   variable = disp_z
  #   use_displaced_mesh = false
  #   component = 1
  # []

  [poro_r]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_r
    use_displaced_mesh = false
    component = 0
  []
  # [poro_z]
  #   type = PorousFlowEffectiveStressCoupling
  #   variable = disp_z
  #   use_displaced_mesh = false
  #   component = 1
  # []
[]

[AuxKernels]
  [rate]
    type = FunctionAux
    variable = rate
    execute_on = timestep_begin
    function = decay_rate
  []
  [pgas]
    type = PorousFlowPropertyAux
    property = pressure
    phase = 1
    variable = pgas
  []
  [sLi]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 0
    variable = sLi
  []
  [stress_rr]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_rr
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
[]

[Functions]
  [decay_rate]
    # Eqn(26) of the first paper of LaForce et al.
    # Ka * (rho C)_a = 10056886.914
    # h = 11
    type = ParsedFunction
    value = 'sqrt(10056886.914/t)/11.0'
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'temp pli gas_saturation disp_r'
    number_fluid_phases = 2
    number_fluid_components = 2
  []
  [pc]
    type = PorousFlowCapillaryPressureConst
    pc = 0
  []
[]

[Modules]
  [FluidProperties]
    [li]
      type = SimpleFluidProperties
      bulk_modulus = 11E9
      viscosity = 1.0e-2
      density0 = 534
      cp = 4194
      cv = 4186
      porepressure_coefficient = 0.0
      thermal_expansion = 0
    []
    [gas]
      type = SimpleFluidProperties
      bulk_modulus = 11e9
      density0 = 516.48
      viscosity = 0.001e-3
      cv = 2920.5
      cp = 2920.5
      porepressure_coefficient = 0.0
      thermal_expansion = 0
    []
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = temp
  []
  [ppss]
    type = PorousFlow2PhasePS
    phase0_porepressure = pli
    phase1_saturation = gas_saturation
    capillary_pressure = pc
  []
  [massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = 'massfrac_ph0_sp0 massfrac_ph1_sp0'
  []
  [water]
    type = PorousFlowSingleComponentFluid
    fp = li
    phase = 0
  []
  [gas]
    type = PorousFlowSingleComponentFluid
    fp = gas
    phase = 1
  []
  [porosity_reservoir]
    type = PorousFlowPorosityConst
    porosity = 0.2
  []
  [permeability_reservoir]
    type = PorousFlowPermeabilityConst
    permeability = '1e-12 0 0  0 0 0  0 0 0'
  []
  [relperm_liquid]
    type = PorousFlowRelativePermeabilityCorey
    n = 4
    phase = 0
    s_res = 0.200
    sum_s_res = 0.405
  []
  [relperm_gas]
    type = PorousFlowRelativePermeabilityBC
    phase = 1
    s_res = 0.205
    sum_s_res = 0.405
    nw_phase = true
    lambda = 2
  []
  [thermal_conductivity_reservoir]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '0 0 0  0 1.320 0  0 0 0'
    wet_thermal_conductivity = '0 0 0  0 1.320 0  0 0 0'
  []
  [internal_energy_reservoir]
    type = PorousFlowMatrixInternalEnergy
    specific_heat_capacity = 1100
    density = 2350.0
  []
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    shear_modulus = 5.0E9
    poissons_ratio = 0.25
  []
  [strain]
    type = ComputeAxisymmetricRZSmallStrain
    eigenstrain_names = 'thermal_contribution'
  []
  # [ini_strain]
  #   type = ComputeEigenstrainFromInitialStress
  #   initial_stress = '-12.8E6 0 0  0 -51.3E6 0  0 0 -12.8E6'
  #   eigenstrain_name = ini_stress
  # []
  [thermal_contribution]
    type = ComputeThermalExpansionEigenstrain
    temperature = temp
    stress_free_temperature = 298
    thermal_expansion_coeff = 0.0
    eigenstrain_name = thermal_contribution
  []
  [stress]
    type = ComputeLinearElasticStress
  []
  [eff_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
  []
  [vol_strain]
    type = PorousFlowVolumetricStrain
  []
[]

[BCs]
  [outer_pressure_fixed]
    type = DirichletBC
    boundary = right
    value = 1.e6
    variable = pli
  []
  [outer_saturation_fixed]
    type = DirichletBC
    boundary = right
    value = 0.0
    variable = gas_saturation
  []
  [outer_temp_fixed]
    type = DirichletBC
    boundary = right
    value = 293
    variable = temp
  []
  [fixed_outer_r]
    type = DirichletBC
    variable = disp_r
    value = 0
    boundary = right
  []
  [li_injection]
    type = PorousFlowSink
    boundary = left
    variable = pli
    use_mobility = false
    use_relperm = false
    fluid_phase = 0
    flux_function = 'if (t <= 36000,-66.413e-9*6.8, 664.13e-9*6.8)'
    # flux_function = 'min(t/100.0,1)*(-2.294001475)' # 5.0E5 T/year = 15.855 kg/s, over area of 2Pi*0.1*11
  []
  [cold_co2]
    type = DirichletBC
    boundary = left
    variable = temp
    value = 293
  []
  [cavity_pressure_x]
    type = Pressure
    boundary = left
    variable = disp_r
    component = 0
    postprocessor = p_bh # note, this lags
    use_displaced_mesh = false
  []
  [co2_injection]
    type = PorousFlowSink
    boundary = left
    variable = gas_saturation
    use_mobility = false
    use_relperm = false
    fluid_phase = 1
    # flux_function = 'min(t/100.0,1)*(-2.294001475e-3)' # 5.0E5 T/year = 15.855 kg/s, over area of 2Pi*0.1*11
    flux_function = '-1e-9'
  []
[]

[Postprocessors]
  [p_bh]
    type = PointValue
    variable = pli
    point = '0.1e-6 0 0'
    execute_on = timestep_begin
    use_displaced_mesh = false
  []
[]

[VectorPostprocessors]
  [ptsuss]
    type = LineValueSampler
    use_displaced_mesh = false
    start_point = '0.12e-6 0 0'
    end_point = '10e-6 0 0'
    sort_by = x
    num_points = 500
    outputs = csv
    variable = 'pli temp gas_saturation disp_r stress_rr stress_tt'
  []
[]

[Preconditioning]
  active = preferred_but_might_not_be_installed
  [smp]
    type = SMP
    full = true
    #petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap '
                          '-snes_atol -snes_rtol -snes_max_it'
    petsc_options_value = 'gmres      asm      lu           NONZERO                   2              '
                          ' 1E2       1E-5        500'
  []
  [mumps]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -pc_factor_shift_type -snes_rtol '
                          '-snes_atol -snes_max_it'
    petsc_options_value = 'lu       mumps                         NONZERO               1E-5       '
                          '1E2       50'
  []
  [preferred_but_might_not_be_installed]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  []
[]

[Executioner]
  type = Transient
  automatic_scaling = true
  solve_type = NEWTON
  end_time = 3600
  resid_vs_jac_scaling_param = 0.5
  # nl_abs_tol = 1e-11
  nl_rel_tol = 1e-06
  #dtmax = 1e6
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 10
    growth_factor = 1.1
  []
[]

[Outputs]
  print_linear_residuals = false
  # sync_times = '3600'# 86400 2.592E6 1.5768E8'
  perf_graph = true
  exodus = true
  [csv]
    type = CSV
    sync_only = false
  []
[]
