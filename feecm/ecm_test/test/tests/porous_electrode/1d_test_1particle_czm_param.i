# Each particle type is in parallel and there is no individual contact resistance
# So the total cell potential is the same
# We should be able to monitor utilization
# all the cell voltage variables can be combined into a single variable
# ---- Optimal Units for Electro-chemo-mechanical Calculations ------ #

# -------------- Base Units ------------------------- #
# Length              = 1e-6 m -> mum
# Charge              = 1e-9 C -> nC
# Amt of substance    = 1e-12 mol -> pmol
# Temperature         = K
# Mass                = kg
# Time                = s
# Voltage             = 1e-3 V -> mV
#
# # -------------- Derived Quantities -------------------#
# Force                = 1e-6 N -> muN
# Stress/Pressure     = 1e6 N/m^2 => MPa
# Energy              = 1e-12 J -> pJ
#
# Concentration       = 1e6 mol/m^3  -> pmol/mum^3
# Concentration Flux  = mol/m^2/s  -> pmol/mum^2/s
# Diffusivity         = 1e-12 m^2/s -> mum^2/s
# Chemical Potential  = J/mol -> pJ/pmol
# Molar Volume        = 1e-6 m^3/mol -> mum^3/pmol
#
#
# Current             = 1e-9 A => nA => nC/s
# Current density     = 1e3 A/m^2 => nA/mum^2
# Resistance          = 1e6 ohm => kgm^2s^-3A^-2 => kg mu^2 s^-3 nA^-2 => Kohm
# Resistivity         = ohm-m
# Electric Field      = 1e3 V/m
# Conductivity        = S/m
# Charge transfer res = 1e-6 ohm m^2
#
# Thermal Conduc      = 1e-6 W/m/K
#
# # ----------------------------------------------------
# Gas Constant        = 8.314 J/mol/K
# Faradays constant   = 96.4853329 nC/pmol
gas_constant = 8.314 # J/mol/K
F_const = 96.4853329 # nC/pmol
temperature = 298 # K
solid_phase_conductivity = 100 # S/m
L_anode = 20 #mum
L_electrolyte = 70 #mum
L_cathode = 120 #mum
particle1_cs_init = 11333.32 #mol/m^3
nmc_c_max = 50060 # mol/m^3
nmc_c_min = 11333.32 # mol/m^3
nmc_exchange_current_density = 5.43e-2 #A/m^2
particle1_size = 15.0 # um
total_volume_fraction = 0.6
particle1_volume_fraction = 0.6
cutoff_voltage = 3.1 #V
cathode_electronic_conductivity = 200.0 #S/m
anode_electronic_conducitivity = 200.0 # S/m
cathode_ionic_conducitivity = 1 #mS/cm
anode_ionic_conductivity = 1 #mS/cm
cathode_solid_phase_diffusivity = 1e-14 #m^2/s
cathode_potential_file = "nmc_2_eq_pot.csv"
cathode_capacity = 200 #mAh/g
n_elements = 210
elem_type = "EDGE3"
cathode_socmin = 0.222
cathode_socmax = 0.942
C_rate = 1.0

[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = 0
    xmax = '${fparse L_anode + L_electrolyte + L_cathode}'
    nx = ${n_elements}
    elem_type = ${elem_type}
    # ny = 2
  []
  [anode]
    type = ParsedSubdomainMeshGenerator
    input = msh
    combinatorial_geometry = 'x <= ${L_anode}'
    block_id = 100
    block_name = 'Anode'
  []
  [cathode]
    type = ParsedSubdomainMeshGenerator
    input = anode
    combinatorial_geometry = 'x >= ${fparse L_anode + L_electrolyte}'
    block_id = 300
    block_name = 'Cathode'
  []
  [electrolyte]
    type = ParsedSubdomainMeshGenerator
    input = cathode
    block_id = 200
    combinatorial_geometry = 'x > ${L_anode} & x < ${fparse L_anode + L_electrolyte}'
    block_name = 'Electrolyte'
  []
  [break]
    type = BreakMeshByBlockGenerator
    input = electrolyte
    split_interface = true
    show_info = true
    # surrounding_blocks = '200 100'
    block_pairs = '100 200'
  []
[]

[Problem]
  type = ReferenceResidualProblem
  # solution_variables = 'ux uy normal_lm thermal_lm li_ion_V thermal_lm2 li_metal_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  # group_variables = 'phi_s ph_e particle1_ave_c'
  # group_variables = 'ux uy; li_ion_V li:_metal_conc thermal_lm'
  group_variables = 'particle1_ave_c particle1_cs'
  acceptable_iterations = 2
  # restart_file_base = check/full_model_cp/0002
[]

[Variables]
  [phi_e]
    # order = FIRST
    block = 'Cathode Electrolyte'
  []
  [phi_s]
    block = 'Cathode Anode'
    # order = FIRST
    # initial_condition = 4000.0
  []
  [particle1_ave_c]
    block = 'Cathode'
    # order = SECOND
    initial_condition = '${units ${particle1_cs_init} mol/m^3 -> pmol/mum^3}'
  []
  [particle1_cs]
    block = 'Cathode'
    initial_condition = '${units ${particle1_cs_init} mol/m^3 -> pmol/mum^3}'
  []
[]

[ICs]
  [V_Cathode]
    type = ConstantIC
    variable = phi_s
    value = '${units 4 V -> mV}'
    block = Cathode
  []
  [V_Anode]
    type = ConstantIC
    variable = phi_s
    value = 0.0
    block = Anode
  []
[]

[AuxVariables]
  [ionic_current_x]
    order = FIRST
    family = MONOMIAL
  []

  [ionic_current_y]
    order = FIRST
    family = MONOMIAL
  []

  [electronic_current_x]
    order = FIRST
    family = MONOMIAL
    block = 'Cathode Anode'
  []
  [electronic_current_y]
    order = FIRST
    family = MONOMIAL
    block = 'Cathode Anode'
  []
[]

[AuxKernels]

  [i_s_x]
    type = DiffusionFluxAux
    variable = electronic_current_x
    diffusion_variable = phi_s
    diffusivity = electronic_conductivity
    component = x
    block = 'Cathode Anode'
  []
  [i_s_y]
    type = DiffusionFluxAux
    variable = electronic_current_y
    diffusion_variable = phi_s
    diffusivity = electronic_conductivity
    component = y
    block = 'Cathode Anode'
  []
  [i_e_x]
    type = DiffusionFluxAux
    variable = ionic_current_x
    diffusion_variable = phi_e
    diffusivity = electrolyte_conductivity
    component = x
    block = 'Cathode Electrolyte'
  []
  [i_e_y]
    type = DiffusionFluxAux
    variable = ionic_current_y
    diffusion_variable = phi_e
    diffusivity = electrolyte_conductivity
    component = y
    block = 'Cathode Electrolyte'
  []
[]

[Functions]
  [eq_pot]
    type = PiecewiseLinear
    # data_file = 'nmc_811_eq_pot.csv'
    data_file = ${cathode_potential_file}
    format = columns
    scale_factor = '${units 1 V -> mV}'
    extrap = true
  []
  [cathode_areal_capacity]
    type = ParsedFunction
    value = 'ceqref * (socmax-socmin)* vol_fraction * L * F / 3600.0'
    vars = 'ceqref socmax socmin vol_fraction L F'
    vals = '${units ${nmc_c_max} mol/m^3 -> pmol/mum^3} ${cathode_socmax} ${cathode_socmin} '
           '${total_volume_fraction} ${L_cathode} ${F_const}'
  []
  [current_density]
    type = ParsedFunction
    value = '1.0 * cap * crate'
    vars = 'cap crate'
    vals = 'cathode_areal_capacity ${C_rate}'
  []
[]

[Kernels]
  [particle1_surface_concentration]
    type = SurfaceConcentration2ParameterAppox
    base_name = 'particle1'
    variable = particle1_cs
    average_concentration_var = particle1_ave_c
    flux_property = butler_volmer_flux
    solid_diffusivity = D_s
    particle_size = particle_size
    block = 'Cathode'
    # scale = 0.2
    # extra_vector_tags = 'ref'
  []

  [V_electrolyte]
    type = ADMatDiffusion
    variable = phi_e
    diffusivity = 'electrolyte_conductivity'
    # thermal_conductivity = thermal_conducitivity
  []
  [V_particle1]
    type = ADMatDiffusion
    variable = phi_s
    block = 'Cathode Anode'
    diffusivity = 'electronic_conductivity'
  []

  [source_electrolyte]
    type = ADMultipleButlerVolmerForce
    variable = phi_e
    scale = -1.0
    materials = 'particle1_bulter_volmer_material'
    mat_prop_name = "butler_volmer_current_force"
    block = 'Cathode'
    extra_vector_tags = 'ref'
  []

  [source_phi]
    type = ADMultipleButlerVolmerForce
    variable = phi_s
    scale = -1.0
    materials = 'particle1_bulter_volmer_material'
    mat_prop_name = "butler_volmer_current_force"
    block = 'Cathode'
    extra_vector_tags = 'ref'
  []

  [source_particle1_ave_c]
    type = ADButlerVolmerForce
    base_name = 'particle1'
    variable = particle1_ave_c
    mat_prop_name = "butler_volmer_flux_force"
    block = 'Cathode'
    scale = '1.0'
    extra_vector_tags = 'ref'
  []

  [cs_time_parcticle1]
    type = ADTimeDerivative
    variable = particle1_ave_c
    block = 'Cathode'
    # extra_vector_tags = 'ref'
  []
[]

[InterfaceKernels]
  [interface_anode_electrolyte]
    type = CZMDiffusiveVariableKernelSmallStrain
    boundary = 'Anode_Electrolyte'
    variable = phi_s
    neighbor_var = phi_e
    include_gap = false
    extra_vector_tags = 'ref'
    flux_global_name = flux_global
    component = 0
  []
[]

[Materials]
  [interface_material]
    type = CZMButlerVolmer
    boundary = 'Anode_Electrolyte'
    surfaceType = SECONDARY
    conductanceType = EXCHANGE_CURRENT_DENSITY
    computeType = BUTLER_VOLMER
    max_separation = 100
    k_function = '1e-1'
    variable = phi_s
    include_equil = false
    include_gap = false
    temperature = ${temperature}
  []
  [interface_jump]
    type = CZMComputeVariableJumpSmallStrain
    variable = phi_s
    neighbor_variable = phi_e
    boundary = 'Anode_Electrolyte'
    include_gap = false
  []
  [interface_flux]
    type = CZMComputeGlobalFluxSmallStrain
    boundary = 'Anode_Electrolyte'
    include_gap = false
  []

  [ionic_conduction_electrolyte]
    type = ADIsotropicDiffusionMaterial
    diffusion_coef = '${units ${cathode_ionic_conducitivity} mS/cm -> S/m}'
    diffusivity_name = 'electrolyte_conductivity'
    block = 'Electrolyte'
  []
  [ionic_conduction_Cathode]
    type = ADIsotropicDiffusionMaterial
    diffusion_coef = '${units ${cathode_ionic_conducitivity} mS/cm -> S/m}'
    diffusivity_name = 'electrolyte_conductivity'
    block = 'Cathode'
  []
  # [ionic_conduction_Cathode]
  #   type = PorousIsotropicDiffusionMaterial
  #   volume_fraction = 0.15
  #   diffusion_coef = 0.2
  #   Bruggeman_factor = 1.5
  #   diffusion_model = BRUGGEMAN
  #   diffusivity_name = 'electrolyte_conductivity'
  #   block = 'Cathode'
  # []
  [electronic_conduction]
    type = ADIsotropicDiffusionMaterial
    # volume_fraction = 0.85
    diffusion_coef = '${units ${cathode_electronic_conductivity} S/m -> S/m}'
    diffusivity_name = 'electronic_conductivity'
    block = 'Cathode'
  []
  [electronic_conduction_Anode]
    type = ADIsotropicDiffusionMaterial
    diffusion_coef = '${units ${anode_electronic_conducitivity} S/m -> S/m}'
    block = 'Anode'
    diffusivity_name = 'electronic_conductivity'
  []
  [solid_phase_diff]
    type = ADIsotropicDiffusionMaterial
    diffusion_coef = '${units ${cathode_solid_phase_diffusivity} m^2/s -> mum^2/s}'
    diffusivity_name = 'D_s'
    block = 'Cathode'
  []
  [particle1_bulter_volmer_material]
    type = ADButlerVolmerMaterial
    base_name = 'particle1'
    electrode_potential_var = phi_s
    electrolyte_potential_var = phi_e
    gas_constant = ${gas_constant}
    faraday = ${F_const}
    temperature = ${temperature}
    surface_to_volume = surface_to_volume
    block = 'Cathode'
    include_equilibrium = true
    exchange_current_density = exchange_current_density
  []

  [particle1_exchange_current_density]
    type = ADExchangeCurrentDensityMaterial
    base_name = 'particle1'
    exchange_current_density_type = "LI_INSERTION"
    reference_current_density = '${units ${nmc_exchange_current_density} A/m^2 -> nA/mum^2}'
    # reference_current_density = 5.43e-5
    concentration = particle1_cs
    # cmax = '${units 50060 mol/m^3 -> umol / um^3}'
    cmax = '${units ${nmc_c_max} mol/m^3 -> pmol/mum^3}'
    block = 'Cathode'
  []

  [particle1_equilibirum_potential]
    type = ADComputeEquilibriumPotential
    base_name = 'particle1'
    include_reaction_rate = true
    reaction_rate_function = eq_pot
    # cref = '${units 50060 mol/m^3 -> umol/um^3}'
    cref = '${units ${nmc_c_max} mol/m^3 -> pmol/mum^3}'
    concentration = particle1_cs
    block = 'Cathode'
  []

  [particle1_Cathode_particle_size]
    type = ADParticleSize
    base_name = 'particle1'
    particleType = SPHERICAL
    # particle_size = '${units 15.0 mum}'
    particle_size = ${particle1_size}
    volume_fraction = ${total_volume_fraction}
    block = Cathode
  []
[]

[BCs]
  [current]
    type = FunctionNeumannBC
    variable = phi_s
    # function = 'if (t >0 ,-1e-3, 0.0)'
    # function = '-0.69
    # value = '${units 81.145 A/m^2 -> nA/mum^2}'
    function = current_density
    # value = -69.55e-2
    # function = 'if (t < 1, 1e-3*t, 1e-3)'
    boundary = right
    extra_vector_tags = 'ref'
  []

  [V0]
    type = ADDirichletBC
    variable = phi_s
    value = 0
    boundary = left
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]
# [Dampers]
#   [const_damp]
#     type = ConstantDamper
#     damping = 0.5
#     min_damping = 0.1
#   []
# []
[UserObjects]
  [voltage_cutoff1]
    type = Terminator
    expression = 'cell_voltage < ${units ${cutoff_voltage} V -> mV}'
  []

  # [max_conc1]
  #   type = Terminator
  #   expression = 'particle1_avg_soc > 0.972'
  # []

[]
[Postprocessors]
  [particle1_max_soc]
    type = ADElementExtremeMaterialProperty
    mat_prop = "particle1_state_of_charge"
    block = Cathode
    value_type = max
  []

  [particle1_avg_soc]
    type = ADElementAverageMaterialProperty
    block = Cathode
    mat_prop = particle1_state_of_charge
  []

  [cell_voltage]
    type = ElementAverageValue
    variable = phi_s
    block = Cathode
  []

  [Electrolyte_potential]
    type = SideAverageValue
    variable = phi_e
    boundary = 'right'
  []
  [particle1_concentration]
    type = ElementAverageValue
    variable = particle1_cs
    block = Cathode
  []

  [nonlin_it]
    type = NumNonlinearIterations
  []
  [cumulative_nonlin_it]
    type = CumulativeValuePostprocessor
    postprocessor = nonlin_it
  []
  [volume_fraction_integral]
    type = ADElementIntegralMaterialProperty
    mat_prop = "particle1_volume_fraction"
    block = 'Cathode'
  []
  [surface_to_volume_integral]
    type = ADElementIntegralMaterialProperty
    mat_prop = "particle1_surface_to_volume"
    block = 'Cathode'
  []
  [empty]
    type = EmptyPostprocessor
  []
  [volume]
    type = VolumePostprocessor
    block = 'Cathode'
  []
  [cell_capacity]
    type = ParsedPostprocessor
    function = 'empty + t * ${cathode_capacity} / 3600'
    pp_names = 'empty'
    use_t = true
  []
  [applied_current_density]
    type = FunctionValuePostprocessor
    function = current_density
  []
[]

[Executioner]
  # Executioner
  type = Transient
  # normalize_solution_diff_norm_by_dt = false
  solve_type = 'NEWTON'
  line_search = basic
  petsc_options = '-snes_converged_reason -ksp_converged_reason -pc_svd_monitor '
                  '-snes_linesearch_monitor'

  petsc_options_iname = '-pc_type -pc_mat_solver_package -pc_factor_shift_type '
                        '-pc_factor_shift_amount' # -mat_mffd_err'
  petsc_options_value = 'lu       mumps NONZERO               1e-15                   ' #1e-10'
  # petsc_options_iname = '-pc_type -pc_mat_solver_package -pc_factor_shift_type -pc_factor_shift_amount'
  # petsc_options_value = 'lu mumps NONZERO               1e-15                   1e-5'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  l_max_its = 20
  start_time = 0.0
  dt = 0.1
  dtmin = 1e-6
  dtmax = 200
  # end_time = 360000
  end_time = 1
  nl_max_its = 25
  num_steps = 5
  automatic_scaling = true
  resid_vs_jac_scaling_param = 0.5
  compute_scaling_once = false
  # scaling_group_variables = 'phi_s phi_e particle1_ave_c'
  # off_diagonals_in_auto_scaling = true
  # scheme = 'BDF2'
  # [TimeStepper]
  #   type = AB2PredictorCorrector
  #   dt = 0.005
  #   max_increase = 2.0
  #   e_max = 10
  #   e_tol = 1
  # []
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.005
    growth_factor = 2.0
    cutback_factor = 0.5
    optimal_iterations = 25
  []
  # [TimeIntegrator]
  #   type = BDF2
  # []
[]

[Outputs]
  # file_base = '1d_2_4'
  exodus = true
  csv = true
  [ex]
    type = Exodus
    # file_base = '1d_2_4_ex'
    output_material_properties = true
    # execute_on = 'INITIAL LINEAR NONLINEAR'
  []
[]
[Controls]
  [start_c]
    type = TimePeriod
    disable_objects = 'Kernel::cs_time_parcticle1 Kernel::source_particle1_ave_c '
                      'Kernel::particle1_surface_concentration'
    start_time = '0'
    end_time = '0.005'
  []
[]
