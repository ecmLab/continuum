[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = 0
    xmax = 50
    nx = 30
    elem_type = EDGE3
    # ny = 2
  []
  [rename]
    type = RenameBlockGenerator
    input = msh
    old_block = 0
    new_block = 'Electrolyte'
  []
  [subdomain_id]
    type = ParsedSubdomainMeshGenerator
    input = rename
    combinatorial_geometry = 'x > 30'
    block_id = 100
    block_name = 'Cathode'
  []
  [cc]
    type = ParsedSubdomainMeshGenerator
    input = subdomain_id
    combinatorial_geometry = 'x < 1.0'
    block_id = 200
    block_name = 'CC'
  []
[]

[Problem]
  type = ReferenceResidualProblem
  # solution_variables = 'ux uy normal_lm thermal_lm li_ion_V thermal_lm2 li_metal_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  # group_variables = 'phi_s ph_e particle1_ave_c'
  # group_variables = 'ux uy; li_ion_V li:_metal_conc thermal_lm'
  group_variables = 'particle1_ave_c particle1_cs;  particle2_ave_c particle2_cs'
  acceptable_iterations = 2
  # restart_file_base = check/full_model_cp/0002
[]

[Variables]
  [phi_e]
    # order = FIRST
  []
  [particle1_phi]
    block = 'Cathode'
    # order = FIRST
    initial_condition = 4000.0
  []
  [particle2_phi]
    block = 'Cathode'
    # order = FIRST
    initial_condition = 4000.0
  []
  [particle1_ave_c]
    block = 'Cathode'
    # order = SECOND
    initial_condition = 0.1
  []
  [particle1_cs]
    block = 'Cathode'
    initial_condition = 0.1
  []
  [particle2_ave_c]
    block = 'Cathode'
    # order = SECOND
    initial_condition = 0.1
  []
  [particle2_cs]
    block = 'Cathode'
    initial_condition = 0.1
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
    block = 'Cathode'
  []
  [electronic_current_y]
    order = FIRST
    family = MONOMIAL
    block = 'Cathode'
  []
[]

[AuxKernels]

  [i_s_x]
    type = DiffusionFluxAux
    variable = electronic_current_x
    diffusion_variable = particle1_phi
    diffusivity = electronic_conductivity
    component = x
    block = 'Cathode'
  []
  [i_s_y]
    type = DiffusionFluxAux
    variable = electronic_current_y
    diffusion_variable = particle1_phi
    diffusivity = electronic_conductivity
    component = y
    block = 'Cathode'
  []
  [i_e_x]
    type = DiffusionFluxAux
    variable = ionic_current_x
    diffusion_variable = phi_e
    diffusivity = electrolyte_conductivity
    component = x
  []
  [i_e_y]
    type = DiffusionFluxAux
    variable = ionic_current_y
    diffusion_variable = phi_e
    diffusivity = electrolyte_conductivity
    component = y
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
  [particle1_surface_concentration]
    type = SurfaceConcentration2ParameterAppox
    base_name = 'particle1'
    variable = particle1_cs
    average_concentration_var = particle1_ave_c
    flux_property = butler_volmer_flux
    solid_diffusivity = D_s
    particle_size = particle_size
    block = 'Cathode'
    scale = 0.5
    # extra_vector_tags = 'ref'
  []
  [particle2_surface_concentration2]
    type = SurfaceConcentration2ParameterAppox
    base_name = 'particle2'
    variable = particle2_cs
    average_concentration_var = particle2_ave_c
    flux_property = butler_volmer_flux
    solid_diffusivity = D_s
    particle_size = particle_size
    block = 'Cathode'
    scale = 0.5
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
    variable = particle1_phi
    block = 'Cathode'
    diffusivity = 'electronic_conductivity'
  []
  [V_particle2]
    type = ADMatDiffusion
    variable = particle2_phi
    block = 'Cathode'
    diffusivity = 'electronic_conductivity'
  []

  [source_electrolyte]
    type = ADMultipleButlerVolmerForce
    variable = phi_e
    scale = -1.0
    materials = 'particle1_bulter_volmer_material particle2_bulter_volmer_material'
    mat_prop_name = "butler_volmer_current_force"
    block = 'Cathode'
    extra_vector_tags = 'ref'
  []

  [source_particle1_phi]
    type = ADButlerVolmerForce
    base_name = 'particle1'
    variable = particle1_phi
    scale = 1.0
    mat_prop_name = "butler_volmer_current_force"
    block = 'Cathode'
    extra_vector_tags = 'ref'
  []

  [source_particle2_phi]
    type = ADButlerVolmerForce
    base_name = 'particle2'
    variable = particle2_phi
    scale = 1.0
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

  [source_particle2_ave_c]
    type = ADButlerVolmerForce
    base_name = 'particle2'
    variable = particle2_ave_c
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
  [cs_time_parcticle2]
    type = ADTimeDerivative
    variable = particle2_ave_c
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
    block = 'Cathode'
  []
  [solid_phase_diff]
    type = ADIsotropicDiffusionMaterial
    diffusion_coef = 0.05
    diffusivity_name = 'D_s'
    block = 'Cathode'
  []
  [particle1_bulter_volmer_material]
    type = ADButlerVolmerMaterial
    base_name = 'particle1'
    electrode_potential_var = particle1_phi
    electrolyte_potential_var = phi_e
    gas_constant = 8.3145
    faraday = 96.4852239
    temperature = 298
    surface_to_volume = surface_to_volume
    block = 'Cathode'
    include_equilibrium = true
    exchange_current_density = exchange_current_density
  []

  [particle2_bulter_volmer_material]
    type = ADButlerVolmerMaterial
    base_name = 'particle2'
    electrode_potential_var = particle2_phi
    electrolyte_potential_var = phi_e
    gas_constant = 8.3145
    faraday = 96.4852239
    temperature = 298
    surface_to_volume = surface_to_volume
    block = 'Cathode'
    include_equilibrium = true
    exchange_current_density = exchange_current_density
  []

  [particle1_exchange_current_density]
    type = ADExchangeCurrentDensityMaterial
    base_name = 'particle1'
    exchange_current_density_type = "LI_INSERTION"
    reference_current_density = 4.5e-2
    concentration = particle1_cs
    cmax = 0.50060
    block = 'Cathode'
  []

  [particle2_exchange_current_density]
    type = ADExchangeCurrentDensityMaterial
    base_name = 'particle2'
    exchange_current_density_type = "LI_INSERTION"
    reference_current_density = 4.5e-2
    concentration = particle2_cs
    cmax = 0.50060
    block = 'Cathode'
  []

  [particle1_equilibirum_potential]
    type = ADComputeEquilibriumPotential
    base_name = 'particle1'
    include_reaction_rate = true
    reaction_rate_function = eq_pot
    cref = 0.50060
    concentration = particle1_cs
    block = 'Cathode'
  []
  [particle2_equilibirum_potential]
    type = ADComputeEquilibriumPotential
    base_name = 'particle2'
    include_reaction_rate = true
    reaction_rate_function = eq_pot
    cref = 0.50060
    concentration = particle2_cs
    block = 'Cathode'
  []

  [particle1_Cathode_particle_size]
    type = ADParticleSize
    base_name = 'particle1'
    particle_size = 0.25
    volume_fraction = 0.5
    block = Cathode
  []
  [particle2_Cathode_particle_size]
    type = ADParticleSize
    base_name = 'particle2'
    particle_size = 0.25
    volume_fraction = 0.5
    block = Cathode
  []

[]

[BCs]
  [particle1_current]
    type = FunctionNeumannBC
    variable = particle1_phi
    # function = 'if (t >0 ,-1e-3, 0.0)'
    function = '-1.2e-1'
    # function = 'if (t < 1, 1e-3*t, 1e-3)'
    boundary = right
    extra_vector_tags = 'ref'
  []
  [particle2_current]
    type = FunctionNeumannBC
    variable = particle2_phi
    # function = 'if (t >0 ,-1e-3, 0.0)'
    function = '-1.2e-1'
    # function = 'if (t < 1, 1e-3*t, 1e-3)'
    boundary = right
    extra_vector_tags = 'ref'
  []

  # [Current]
  #   type = FunctionNeumannBC
  #   variable = phi_e
  #   function = '2.4e-1'
  #   boundary = left
  #   extra_vector_tags = 'ref'
  # []
  [V0]
    type = ADDirichletBC
    variable = phi_e
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
    expression = 'particle1_voltage < 3000'
  []
  [voltage_cutoff2]
    type = Terminator
    expression = 'particle2_voltage < 3000'
  []

  [max_conc1]
    type = Terminator
    expression = 'particle1_max_soc > 0.972'
  []
  [max_conc2]
    type = Terminator
    expression = 'particle2_max_soc > 0.972'
  []

[]
[Postprocessors]
  [particle1_max_soc]
    type = ADElementExtremeMaterialProperty
    mat_prop = "particle1_state_of_charge"
    block = Cathode
    value_type = max
  []
  [particle2_max_soc]
    type = ADElementExtremeMaterialProperty
    mat_prop = "particle2_state_of_charge"
    block = Cathode
    value_type = max
  []

  [particle1_voltage]
    type = ElementAverageValue
    variable = particle1_phi
    block = Cathode
  []
  [particle2_voltage]
    type = ElementAverageValue
    variable = particle2_phi
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
  [particle2_concentration]
    type = ElementAverageValue
    variable = particle2_cs
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
  line_search = basic
  petsc_options = '-snes_converged_reason -ksp_converged_reason -pc_svd_monitor '
                '-snes_linesearch_monitor'

  petsc_options_iname = '-pc_type -pc_mat_solver_package -pc_factor_shift_type -pc_factor_shift_amount'# -mat_mffd_err'
  petsc_options_value = 'lu       mumps NONZERO               1e-15                   '#1e-10'
  # petsc_options_iname = '-pc_type -pc_mat_solver_package -pc_factor_shift_type -pc_factor_shift_amount'
  # petsc_options_value = 'lu mumps NONZERO               1e-15                   1e-5'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  l_max_its = 20
  start_time = 0.0
  dt = 0.1
  dtmin = 1e-6
  dtmax = 25
  end_time = 360000
  nl_max_its = 25
  # num_steps = 2
  automatic_scaling = true
  resid_vs_jac_scaling_param = 0.5
  compute_scaling_once = false
  # scaling_group_variables = 'phi_s phi_e particle1_ave_c'
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
  file_base = '1d_2_4'
  exodus = true
  csv = true
  [ex]
    type = Exodus
    file_base = '1d_2_4_ex'
    output_material_properties = true
    # execute_on = 'INITIAL LINEAR NONLINEAR'
  []
[]
