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
    combinatorial_geometry = 'x > 25'
    block_id = 100
    block_name = 'Cathode'
  []
[]

[Problem]
  type = ReferenceResidualProblem
  # solution_variables = 'ux uy normal_lm thermal_lm li_ion_V thermal_lm2 li_metal_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  # group_variables = 'ux uy; li_ion_V li_metal_conc thermal_lm'
  # group_variables = 'disp_'
  acceptable_iterations = 2
  # restart_file_base = check/full_model_cp/0002
[]

[Variables]
  [phi_e]
    order = SECOND
  []
  [phi_s]
    block = Cathode
    order = SECOND
    initial_condition = 4000.0
  []
  [c_s]
    block = Cathode
    order = SECOND
    initial_condition = 0.1
  []
[]

[Functions]
  [eq_pot]
    type = PiecewiseLinear
    data_file = 'nmc_equilibrium_potential.csv'
    format = columns
  []
[]

[Kernels]
  [V]
    type = ADHeatConduction
    variable = phi_e
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
    scale = 1.0
    mat_prop_name = "butler_volmer_current_force"
    block = 'Cathode'
  []
  [source2]
    type = ADCoupledMatForce
    variable = phi_s
    scale = 1.0
    mat_prop_name = "butler_volmer_current_force"
    block = 'Cathode'
  []
  [source3]
    type = ADCoupledMatForce
    variable = c_s
    mat_prop_name = "butler_volmer_flux"
    block = 'Cathode'
    scale = '-0.3'
  []

  [cs_time]
    type = ADTimeDerivative
    variable = c_s
    block = 'Cathode'
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
  [bulter_volmer_material]
    type = ADButlerVolmerConstantExchangeCurrent
    electrode_potential_var = phi_s
    electrolyte_potential_var = phi_e
    gas_constant = 8.3145
    faraday = 96.4852239
    temperature = 298
    particle_size = 10.0
    block = 'Cathode'
    include_equilibrium = true
    exchange_current_density = 1e-3
  []
  [equilibirum_potential]
    type = ADComputeEquilibriumPotential
    include_reaction_rate = true
    reaction_rate_function = eq_pot
    cref = 0.50060
    concentration = c_s
    block = 'Cathode'
  []
[]

[BCs]
  [current]
    type = FunctionNeumannBC
    variable = phi_s
    # function = 'if (t >0 ,-1e-3, 0.0)'
    function = '1e-3'
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

[Executioner]
  # Executioner
  type = Transient

  solve_type = 'NEWTON'
  line_search = none
  petsc_options_iname = '-pc_type -pc_mat_solver_package'
  petsc_options_value = 'lu mumps'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-11
  l_max_its = 20
  start_time = 0.0
  dt = 0.1
  dtmin = 1e-6
  dtmax = 200
  end_time = 3600
  nl_max_its = 200
  num_steps = 100
  automatic_scaling = true
  resid_vs_jac_scaling_param = 0
  compute_scaling_once = false
  # scaling_group_variables = 'disp_x V'
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.005
    growth_factor = 2
    cutback_factor = 0.5
    optimal_iterations = 6
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
    boundary = 'left'
  []
  [Cell_concentration]
    type = ElementAverageValue
    variable = c_s
    block = Cathode
  []
[]

[Outputs]
  exodus = true
  [ex]
    type = Exodus
    output_material_properties = true
    # execute_on = 'INITIAL LINEAR NONLINEAR'
  []
[]
