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
    file = 'test1.msh'
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
# [GlobalParams]
#   displacements = 'ux uy'
# []

[Variables]
  [./V]
  [../]
  [./li_metal_conc]
    initial_condition = 7.35e-3
    block = 'Cathode'
  [../]
  [./thermal_lm]
    block = 'secondary_contact_block'
  [../]
[]
[ICs]
  [./Voltage_cathode]
    type = ConstantIC
    block = 'Cathode'
    value = 3000.0
    variable = V
  [../]
  [./Voltage_electrolyte]
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
  [../]
[]

[Functions]
  [./reaction_rate]
    type = PiecewiseLinear
    data_file = 'nmc_equilibrium_potential.csv'
    format = columns
  [../]
[]

[AuxVariables]
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
  [./bnd_li_flux]
    type = DiffusionFluxNormalToBoundaryAux
    variable = bndliflux
    boundary = 'Electrolyte_bottom Cathode_top'
    diffusion_variable = V
    diffusivity = thermal_conductivity
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
  [./li_ion_flux_z]
    type = ADDiffusionFluxAux
    variable = flux_z
    component = z
    diffusion_variable = V
    diffusivity = thermal_conductivity
    block = 'Cathode Electrolyte'
  [../]
[]


[Kernels]
  [./li_ion_V]
    type = ADHeatConduction
    variable = V
    use_displaced_mesh = false
    block = 'Cathode Electrolyte'
  [../]

  [./li_metal2]
    type = ADMatDiffusion
    variable = li_metal_conc
    diffusivity = diffusivity
    use_displaced_mesh = false
    block = 'Cathode'
  [../]
  [./li_metal_dt]
    type = ADTimeDerivative
    variable = li_metal_conc
    use_displaced_mesh = false
    block = 'Cathode'
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
  [./diffusivity_Li]
    type = ADGenericConstantMaterial
    prop_names = 'diffusivity'
    prop_values = '5'
    block = 'Cathode'
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
    block = Cathode
  [../]

[]

[BCs]
  [./current]
    type = ADNeumannBC
    boundary = 'Cathode_bottom'
    value = -20e-3
    variable = V
    extra_vector_tags = 'ref'
  [../]
  # [./current]
  #   type = ADButlerVolmerBC
  #   current_density = 1e-3
  #   exchange_current_density = 1e-1
  #   faraday = 96.4853
  #   Temperature = 298
  #   R = 8.314462681
  #   variable = V
  #   boundary = 'Cathode_bottom'
  #   # extra_vector_tags = 'ref'
  # [../]
  [./conc]
    type = ScaledCoupledVarNeumannBC
    variable = li_metal_conc
    v = bndliflux
    scale = -1.036428e-2
    # value = 1.03642813e-5
    boundary = 'Cathode_top'
    extra_vector_tags = 'ref'
  [../]
  [./OV]
    type = ADDirichletBC
    boundary = 'Electrolyte_top'
    value = 0.0
    variable = V
  [../]

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
  [./conc]
    type = SideAverageValue
    boundary = 'Cathode_bottom'
    variable = li_metal_conc
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  automatic_scaling = true
  # compute_scaling_once = false
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1     NONZERO               1e-20               '
  dt = 200
  # num_steps = 20
  nl_max_its = 15
  nl_abs_tol = 1e-3
  nl_rel_tol = 1e-6
  dtmax = 1000
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 10
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 15
    # timestep_limiting_postprocessor = matl_ts_min
  [../]
  # end_time = 36000

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
[]
[Outputs]
  exodus = true
  csv = true
[]
