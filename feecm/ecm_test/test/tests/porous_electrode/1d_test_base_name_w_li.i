[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = 0
    xmax = 60
    nx = 40
    elem_type = EDGE3
    # ny = 2
  []
  [anode]
    type = ParsedSubdomainMeshGenerator
    input = msh
    combinatorial_geometry = 'x < 20'
    block_id = 100
    block_name = 'Anode'
  []
  [cathode]
    type = ParsedSubdomainMeshGenerator
    input = anode
    combinatorial_geometry = 'x > 40'
    block_id = 300
    block_name = 'Cathode'
  []
  [electrolyte]
    type = ParsedSubdomainMeshGenerator
    input = cathode
    block_id = 200
    combinatorial_geometry = 'x > 20 & x < 40'
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
  # group_variables = 'phi_s_s ph_e ave_c'
  # group_variables = 'ux uy; li_ion_V li:_metal_conc thermal_lm'
  group_variables = 'ave_c cs'
  acceptable_iterations = 2
  # restart_file_base = check/full_model_cp/0002
[]

[Variables]
  [phi_e]
    # order = FIRST
    # block = 'Cathode Electrolyte'
  []
  [phi_s]
    block = 'Cathode Anode'
    # order = FIRST
    # initial_condition = 4000.0
  []
  [ave_c]
    block = 'Cathode'
    # order = SECOND
    initial_condition = 0.1
  []
  [cs]
    block = 'Cathode'
    initial_condition = 0.1
  []
[]
[ICs]
  [V_Cathode]
    type = ConstantIC
    variable = phi_s
    value = 4000
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
  [disp_x]
  []
[]

[AuxKernels]
  [dummy_x]
    type = ConstantAux
    variable = disp_x
    value = 0
  []
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
    data_file = 'nmc_811_eq_pot.csv'
    format = columns
    scale_factor = 1000.0
  []
[]

[Kernels]
  [surface_concentration]
    type = SurfaceConcentration2ParameterAppox
    #base_name = 'particle1'
    variable = cs
    average_concentration_var = ave_c
    flux_property = butler_volmer_flux
    solid_diffusivity = D_s
    particle_size = particle_size
    block = 'Cathode'
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
    #base_name = 'particle1'
    materials = 'bulter_volmer_material'
    variable = phi_e
    scale = -1.0
    mat_prop_name = "butler_volmer_current_force"
    block = 'Cathode'
    extra_vector_tags = 'ref'
  []

  [source_phi_s]
    type = ADButlerVolmerForce
    #base_name = 'particle1'
    variable = phi_s
    scale = 1.0
    mat_prop_name = "butler_volmer_current_force"
    block = 'Cathode'
    extra_vector_tags = 'ref'
  []

  [source_ave_c]
    type = ADButlerVolmerForce
    #base_name = 'particle1'
    variable = ave_c
    mat_prop_name = "butler_volmer_flux_force"
    block = 'Cathode'
    scale = '-1.0'
    extra_vector_tags = 'ref'
  []

  [cs_time_parcticle1]
    type = ADTimeDerivative
    variable = ave_c
    block = 'Cathode'
    # extra_vector_tags = 'ref'
  []

[]

[Materials]
  [ionic_conduction]
    type = ADIsotropicDiffusionMaterial
    diffusion_coef = 0.02
    diffusivity_name = 'electrolyte_conductivity'
  []
  [electronic_conduction]
    type = ADIsotropicDiffusionMaterial
    diffusion_coef = 1.0
    diffusivity_name = 'electronic_conductivity'
    block = 'Cathode Anode'
  []
  [solid_phase_diff]
    type = ADIsotropicDiffusionMaterial
    diffusion_coef = 0.05
    diffusivity_name = 'D_s'
    block = 'Cathode'
  []
  [bulter_volmer_material]
    type = ADButlerVolmerMaterial
    #base_name = 'particle1'
    electrode_potential_var = phi_s
    electrolyte_potential_var = phi_e
    gas_constant = 8.3145
    faraday = 96.4852239
    temperature = 298
    surface_to_volume = surface_to_volume
    block = 'Cathode'
    include_equilibrium = true
    exchange_current_density = exchange_current_density
  []

  [exchange_current_density]
    type = ADExchangeCurrentDensityMaterial
    #base_name = 'particle1'
    exchange_current_density_type = "LI_INSERTION"
    reference_current_density = 4.5e-2
    concentration = cs
    cmax = 0.50060
    block = 'Cathode'
  []

  [equilibirum_potential]
    type = ADComputeEquilibriumPotential
    #base_name = 'particle1'
    include_reaction_rate = true
    reaction_rate_function = eq_pot
    cref = 0.50060
    concentration = cs
    block = 'Cathode'
  []

  [Cathode_particle_size]
    type = ADParticleSize
    #base_name = 'particle1'
    particle_size = 0.5
    volume_fraction = 1.0
    block = Cathode
  []

[]

[BCs]
  # [Current]
  #   type = FunctionNeumannBC
  #   variable = phi_e
  #   function = '2.4e-1'
  #   boundary = left
  #   extra_vector_tags = 'ref'
  # []
  # [V0]
  #   type = ADDirichletBC
  #   variable = phi_e
  #   value = 0
  #   boundary = right
  # []
  [current]
    type = FunctionNeumannBC
    variable = phi_s
    # function = 'if (t >0 ,-1e-3, 0.0)'
    function = '-2.4e-1'
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
  [voltage_cutoff]
    type = Terminator
    expression = 'Cell_voltage < 3000'
  []
  [max_conc]
    type = Terminator
    expression = 'max_soc > 0.972'
  []
[]
[Postprocessors]
  [max_soc]
    type = ADElementExtremeMaterialProperty
    mat_prop = "particle1_state_of_charge"
    block = Cathode
    value_type = max
  []
  [Cell_voltage]
    type = ElementAverageValue
    variable = phi_s
    block = Cathode
  []
  [Electrolyte_potential]
    type = SideAverageValue
    variable = phi_e
    boundary = 'right'
  []
  [Cell_concentration]
    type = ElementAverageValue
    variable = ave_c
    block = Cathode
  []
  [nonlin_it]
    type = NumNonlinearIterations
  []
  [cumulative_nonlin_it]
    type = CumulativeValuePostprocessor
    postprocessor = nonlin_it
  []
[]

[Executioner]
  # Executioner
  type = Transient
  # normalize_solution_diff_norm_by_dt = false
  solve_type = 'NEWTON'
  line_search = none
  petsc_options = '-snes_converged_reason -ksp_converged_reason -pc_svd_monitor '
                  '-snes_linesearch_monitor'
  petsc_options_iname = '-pc_type -mat_solver_package -pc_factor_shift_type -pc_factor_shift_amount '
  #'-mat_mffd_err'
  petsc_options_value = 'lu   mumps    NONZERO               1e-15                   ' #1e-5'
  nl_rel_tol = 1e-3
  nl_abs_tol = 1e-6
  l_max_its = 20
  start_time = 0.0
  dt = 0.1
  dtmin = 1e-6
  dtmax = 25
  end_time = 360000
  # end_time = 1
  nl_max_its = 25
  # num_steps = 10
  automatic_scaling = true
  resid_vs_jac_scaling_param = 0.5
  compute_scaling_once = false
  # scaling_group_variables = 'phi_s_s phi_e ave_c'
  # off_diagonals_in_auto_scaling = true
  scheme = 'BDF2'
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