# Test for full cell
# Geometry consists of
# 1) Cathode layer with given OCV curve = nmc_equilibrium_potential.csv
# 2) Solid Electrolyte
# 3) Anode at 0 V
# No mechanics included here for simplicity
# The goal is to just check to see of the whole electro-chemistry works
# All the GapDisplacementConductanceConstraint is just Linear Butler Volmer Kinetics
# All voltages are in mV => 3V = 3000 mV
# --- Issues ----
# 1) Changing the charge transfer resistance changes the voltage profile
# 2) Convergence changes when changing R_ct
# Note k in GapDisplacementConductanceConstraint is 1/Rct
[Mesh]
  [./mesh]
    type = FileMeshGenerator
    file = 'test1_quad3.msh'
  [../]
  [./primary_block]
    type = LowerDBlockFromSidesetGenerator
    input = mesh
    sidesets = 'Electrolyte_bottom'
    new_block_name = 'primary_contact_block'
  [../]
  [./secondary_block]
    type = LowerDBlockFromSidesetGenerator
    input = primary_block
    sidesets = 'Cathode_top'
    new_block_name = 'secondary_contact_block'
  [../]

[]

[Variables]
  [./V]
    block = 'Cathode Electrolyte'
  [../]
  [./thermal_lm]
    block = 'secondary_contact_block'
  [../]
  # [./li_metal_conc]
  #   initial_condition = 2.45e-3
  #   block = 'Cathode'
  # [../]
[]
[ICs]
  [./Voltage_cathode]
    type = ConstantIC
    block = 'Cathode'
    value = 3000.0
    variable = V
  [../]
  [./Voltage_elec]
    type = ConstantIC
    block = 'Electrolyte'
    value = 0.0
    variable = V
  [../]

[]

[Constraints]
  [./thermal_constraint]
    type = GapDisplacementConductanceConstraint
    variable = thermal_lm
    secondary_variable = V
    primary_variable = V
    secondary_boundary =  'Cathode_top'
    secondary_subdomain = 'secondary_contact_block'
    primary_boundary = 'Electrolyte_bottom'
    primary_subdomain = 'primary_contact_block'
    # k_function = gapk
    k = 1e-5
    # use_displaced_mesh = true
    displacements = 'ux uy'
    compute_lm_residuals = true
    include_equilibrium_potential = true
    include_gap = false
    extra_vector_tags = 'ref'
  [../]
[]

[Functions]
  [./reaction_rate]
    type = PiecewiseLinear
    data_file = 'nmc_equilibrium_potential.csv'
    format = columns
  [../]
  [./conc]
    type = ParsedFunction
    value = 'if (t >= 100, 7.97252406e-6 * (t-100.0) + 2.45e-3, 2.45e-3)'
  [../]
  [./flux]
    type = ParsedFunction
    value = 'if (t < 100, 0.1e-3*t, 10e-3)'
  [../]
[]

[AuxVariables]
  [./li_metal_conc]
    block = 'Cathode'
  [../]
  [./ux]
    block = 'Cathode Electrolyte'
  [../]
  [./uy]
    block = 'Cathode Electrolyte'
  [../]

  [./flux_x]
    order = FIRST
    family = MONOMIAL
    block = 'Cathode Electrolyte'
  [../]
  [./flux_y]
    order = FIRST
    family = MONOMIAL
    block = 'Cathode Electrolyte'
  [../]
  [./flux_z]
    order = FIRST
    family = MONOMIAL
    block = 'Cathode Electrolyte'
  [../]
  [./bndliflux]
    order = CONSTANT
    family = MONOMIAL
    block = 'Cathode Electrolyte'
  [../]
  [./Eq_pot]
    order = CONSTANT
    family = MONOMIAL
    block = 'Cathode'
  [../]
[]

[AuxKernels]
  [./conc]
    type = FunctionAux
    function = conc
    variable = li_metal_conc
    block = 'Cathode'
  [../]
  [./ux]
    type = ConstantAux
    value = 0
    block = 'Cathode Electrolyte'
    variable = ux
  [../]
  [./uy]
    type = ConstantAux
    value = 0
    block = 'Cathode Electrolyte'
    variable = uy
  [../]

  [./Eq_pot]
    type = ADMaterialRealAux
    variable = Eq_pot
    property = equilibrium_potential
    block = 'Cathode'
  [../]

  [./li_ion_flux_x]
    type = ADDiffusionFluxAux
    variable = flux_x
    component = x
    diffusion_variable = V
    diffusivity = thermal_conductivity
    block = 'Cathode Electrolyte'
  [../]

  [./li_ion_flux_y]
    type = ADDiffusionFluxAux
    variable = flux_y
    component = y
    diffusion_variable = V
    diffusivity = thermal_conductivity
    block = 'Cathode Electrolyte'
  [../]
  [./bnd_li_flux]
    type = DiffusionFluxNormalToBoundaryAux
    variable = bndliflux
    boundary = 'Electrolyte_bottom Cathode_top'
    diffusion_variable = V
    diffusivity = thermal_conductivity
  [../]
[]


[Kernels]
  [./li_ion_V]
    type = ADHeatConduction
    variable = V
    use_displaced_mesh = false
    block = 'Cathode Electrolyte'
  [../]


[]

[Materials]

  [./thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 100
    block = 'Cathode'
  [../]
  [./thermal_conductivity2]
    type = ADHeatConductionMaterial
    thermal_conductivity = 0.1
    block = 'Electrolyte'
  [../]

  [./equilibrium_potential]
    type = ADComputeEquilibriumPotential
    R = 8.31446
    faraday = 96.4853329
    temperature = 298
    cref = 4.9e-2
    concentration = li_metal_conc
    include_conc = false
    include_reaction_rate = true
    # reaction_rate = 780.0
    reaction_rate_function = reaction_rate
    include_mechanical_effects = false
    exclude_elastic_contribution = true
    block = 'Cathode'
  [../]

[]

[BCs]

  [./current]
    type = ADFunctionNeumannBC
    boundary = 'Cathode_bottom'
    function = flux
    variable = V
    extra_vector_tags = 'ref'
  [../]
  [./OV]
    type = ADDirichletBC
    boundary = 'Electrolyte_top'
    value = 0.0
    variable = V
    preset = true
  [../]
  # [./conc]
  #   type = ScaledCoupledVarNeumannBC
  #   variable = li_metal_conc
  #   v = bndliflux
  #   scale = 1.036428e-2
  #   boundary = 'Cathode_top'
  #   extra_vector_tags = 'ref'
  # [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./Voltage_Cathode]
    type = SideAverageValue
    boundary = 'Cathode_bottom'
    variable = V
  [../]

  [./eq_pot]
    type = SideAverageValue
    boundary = 'Cathode_bottom'
    variable = Eq_pot
  [../]
  [./bot_conc]
    type = SideAverageValue
    boundary = 'Cathode_bottom'
    variable = li_metal_conc
  [../]
[]
[UserObjects]
  [./vcutoff]
    type = Terminator
    expression = 'abs(Voltage_Cathode) < 2750'
    execute_on = 'TIMESTEP_END'
  [../]
  [./max_conc]
    type = Terminator
    expression = 'bot_conc > 4.77e-2'
    execute_on = 'TIMESTEP_END'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  automatic_scaling = true
  compute_scaling_once = false
  petsc_options_iname = '-pc_type -pc_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'
  verbose = true
  resid_vs_jac_scaling_param = 0.5
  # petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration -pc_factor_shift_type -pc_factor_shift_amount'
  # petsc_options_value = 'lu superlu_dist basic 1     NONZERO               1e-15               '
  dt = 200
  # line_search = basic
  # num_steps = 20
  nl_max_its = 25
  # nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8
  dtmax = 200
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 10
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 15
    # timestep_limiting_postprocessor = matl_ts_min
  [../]
  # end_time = 10000

  num_steps = 5
  snesmf_reuse_base = true
  scaling_group_variables = 'V thermal_lm'

[]
[Problem]
  type = ReferenceResidualProblem
  # solution_variables = 'ux uy normal_lm thermal_lm li_ion_V thermal_lm2 li_metal_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'V thermal_lm'
  acceptable_iterations = 2
  # coord_type = RZ
[]
[Outputs]
  exodus = true
  csv = true
[]
