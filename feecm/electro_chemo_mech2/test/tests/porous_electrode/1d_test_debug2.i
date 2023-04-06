[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = 0
    xmax = 50
    nx = 50
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
[]

[Problem]
  type = ReferenceResidualProblem
  # solution_variables = 'ux uy normal_lm thermal_lm li_ion_V thermal_lm2 li_metal_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  # group_variables = 'phi_s ph_e c_s'
  # group_variables = 'ux uy; li_ion_V li:_metal_conc thermal_lm'
  # group_variables = 'disp_'
  acceptable_iterations = 2
  # restart_file_base = check/full_model_cp/0002
[]

[Variables]
  [phi_e]
    # order = FIRST
  []
  [phi_s]
    block = Cathode
    # order = FIRST
    initial_condition = 4000.0
  []
  [c_s]
    block = Cathode
    # order = SECOND
    initial_condition = 0.1
  []
[]

[AuxVariables]
  [c_s1]
    order = FIRST
    family = MONOMIAL
    block = 'Cathode'
    # initial_condition = 0.1
  []
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
    block = Cathode
  []
  [electronic_current_y]
    order = FIRST
    family = MONOMIAL
    block = Cathode
  []

  [butler_volmer]
    order = SECOND
    family = MONOMIAL
  []
[]

[AuxKernels]
  [bv]
    type = ADMaterialRealAux
    block = 'Cathode'
    variable = butler_volmer
    property = butler_volmer_current
  []
  [c_s1]
    type = ADMaterialRealAux
    block = 'Cathode'
    variable = c_s1
    property = surface_concentration
  []

  [i_s_x]
    type = DiffusionFluxAux
    variable = electronic_current_x
    diffusion_variable = phi_s
    diffusivity = diffusivity
    component = x
    block = Cathode
  []
  [i_s_y]
    type = DiffusionFluxAux
    variable = electronic_current_y
    diffusion_variable = phi_s
    diffusivity = diffusivity
    component = y
    block = Cathode
  []
  [i_e_x]
    type = DiffusionFluxAux
    variable = ionic_current_x
    diffusion_variable = phi_e
    diffusivity = thermal_conductivity
    component = x
  []
  [i_e_y]
    type = DiffusionFluxAux
    variable = ionic_current_y
    diffusion_variable = phi_e
    diffusivity = thermal_conductivity
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
  [V]
    type = ADHeatConduction
    variable = phi_e
    # thermal_conductivity = thermal_conducitivity
  []
  [V2]
    type = ADMatDiffusion
    variable = phi_s
    block = 'Cathode'
    diffusivity = diffusivity
  []
  [source1]
    type = ADCoupledMatForce
    variable = phi_e
    scale = -1.0
    mat_prop_name = "butler_volmer_current_force"
    block = 'Cathode'
    extra_vector_tags = 'ref'
  []
  [source2]
    type = ADCoupledMatForce
    variable = phi_s
    scale = 1.0
    mat_prop_name = "butler_volmer_current_force"
    block = 'Cathode'
    extra_vector_tags = 'ref'
  []
  [source3]
    type = ADCoupledMatForce
    variable = c_s
    mat_prop_name = "butler_volmer_flux"
    block = 'Cathode'
    scale = '3.0'
    extra_vector_tags = 'ref'
  []

  [cs_time]
    type = ADTimeDerivative
    variable = c_s
    block = 'Cathode'
    # extra_vector_tags = 'ref'
  []

[]

[Materials]
  [ionic_conduction]
    type = ADHeatConductionMaterial
    thermal_conductivity = 0.02
  []
  [electronic_conduction]
    type = ADGenericConstantMaterial
    prop_names = 'diffusivity'
    prop_values = '1.0'
  []
  [solid_phase_diff]
    type = ADGenericConstantMaterial
    prop_names = 'diff_s'
    prop_values = '0.01'
  []
  [bulter_volmer_material]
    # type = ADButlerVolmerConstantExchangeCurrent
    # type = ADButlerVolmerConcentrationDepExchangeCurrent
    type = ADButlerVolmerMaterial
    electrode_potential_var = phi_s
    electrolyte_potential_var = phi_e
    gas_constant = 8.3145
    faraday = 96.4852239
    temperature = 298
    particle_size = 1.0
    block = 'Cathode'
    include_equilibrium = true
    exchange_current_density = exchange_current_density
    # exchange_current_density = 1e-1
  []
  [exchange_current_density]
    type = ADExchangeCurrentDensityMaterial
    exchange_current_density_type = "LI_INSERTION"
    i0 = 4.5e-2
    concentration = c_s
    # surface_concentration = surface_concentration
    cmax = 0.50060
    block = 'Cathode'
  []

  [equilibirum_potential]
    type = ADComputeEquilibriumPotential
    include_reaction_rate = true
    reaction_rate_function = eq_pot
    cref = 0.50060
    concentration = c_s
    # use_surface_concentration = true
    # surface_concentration_prop = surface_concentration
    block = 'Cathode'
  []
  [surface_concentration]
    type = ADComputeSurfaceConcentration
    block = 'Cathode'
    concentration = c_s
    flux = butler_volmer_flux
    particle_size = 1.0
    diffusivity = diff_s
    scale = -1.0
  []
[]

[BCs]
  [current]
    type = FunctionNeumannBC
    variable = phi_s
    # function = 'if (t >0 ,-1e-3, 0.0)'
    function = '-4.8e-1'
    # function = 'if (t < 1, 1e-3*t, 1e-3)'
    boundary = right
    extra_vector_tags = 'ref'
  []
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
  [cutoff]
    type = Terminator
    expression = 'Cell_voltage < 3000'
  []
[]
[Postprocessors]
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
    variable = c_s
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
  petsc_options_iname = '-pc_type -pc_mat_solver_package'
  petsc_options_value = 'lu mumps'
  nl_rel_tol = 1e-3
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
  # scaling_group_variables = 'phi_s phi_e c_s'
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
  file_base = '1d_4_8'
  exodus = true
  csv = true
  [ex]
    type = Exodus
    file_base = '1d_4_8_ex'
    output_material_properties = true
    # execute_on = 'INITIAL LINEAR NONLINEAR'
  []
[]
