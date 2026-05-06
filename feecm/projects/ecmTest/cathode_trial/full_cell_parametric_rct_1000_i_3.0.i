# Test for Si to compare to Bower and Guduru (2014)
# Material Properties
# Molar density of Si = c_0 = 7.874e4 mol/m^3 = 7.874e-2
# Modulus of Si = 100 GPa = 100e3
# poissons_ratio of Si = 0.26
# Vol. expansion of Si (normalized conc ) = 0.7
#    -> Vol. expansion of Si -> 0.7 * 7.874e-2 = omega
# epsilon_0 (strain rate Si) = 0.6e-9 /s
# initial yield stress of Si = 0.12 GPa = 0.12e3
# Stress exponent for plastic flow in Si = 4 -> 0.25
# initial conc of Li in Si = 0.0078 = 0.0078 * 7.87e-2
#applied current density = 0.012 mA/cm^2 = 1.2e-5
# exchange_current_density = 0.001 mA/cm^2 = 1e-6

current_density = 3e-3
k_anode = 1.0e-05
faraday = 96.4853329
temperature = 333
gas_constant = 8.31446
c_ref = 4.9e-2
c_init = 2.45e-3
v_cutoff = 2750
h = 26

[Mesh]
  [./mesh]
    type = FileMeshGenerator
    file = 'test_full4.msh'
  [../]
  [./primary_boundary_cathode]
    type = LowerDBlockFromSidesetGenerator
    input = mesh
    sidesets = 'Electrolyte_bottom'
    new_block_name = 'cathode_p_block'
  [../]
  [./secondary_boundary_cathode]
    type = LowerDBlockFromSidesetGenerator
    input = primary_boundary_cathode
    sidesets = 'Cathode_top'
    new_block_name = 'cathode_s_block'
  [../]
  [./secondary_anode_block]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_boundary_cathode
    sidesets = 'Anode_bottom'
    new_block_name = 'anode_s_block'
  [../]
  [./primary_anode_block]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_anode_block
    sidesets = 'Electrolyte_top'
    new_block_name = 'anode_p_block'
  [../]
[]


[Variables]
  [./V]
  [../]
  [./V_cathode]
    block = 'cathode_s_block'
  [../]
  [./V_anode]
    block = 'anode_s_block'
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
  [./cathode_constraint]
    type = GapDisplacementConductanceConstraint
    variable = V_cathode
    secondary_variable = V
    secondary_boundary =  'Cathode_top'
    secondary_subdomain = 'cathode_s_block'
    primary_boundary = 'Electrolyte_bottom'
    primary_subdomain = 'cathode_p_block'
    # k_function = gapk
    k = 1e-3
    # use_displaced_mesh = true
    displacements = 'ux uy'
    compute_lm_residuals = true
    include_equilibrium_potential = true
  [../]
  [./anode_constraint]
    type = GapDisplacementConductanceConstraint
    variable = V_anode
    secondary_variable = V
    secondary_boundary =  'Anode_bottom'
    secondary_subdomain = 'anode_s_block'
    primary_boundary = 'Electrolyte_top'
    primary_subdomain = 'anode_p_block'
    # k_function = gapk
    k = ${k_anode}
    # use_displaced_mesh = true
    displacements = 'ux uy'
    compute_lm_residuals = true
    include_equilibrium_potential = false
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
    # value = '7.97252406e-6 * t + 2.45e-3'
    value = '${current_density}/${faraday}/${h} * t + ${c_init}'
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
    block = 'Cathode Electrolyte Anode'
  [../]
  [./flux_y]
    order = FIRST
    family = MONOMIAL
    block = 'Cathode Electrolyte Anode'
  [../]
  [./flux_z]
    order = FIRST
    family = MONOMIAL
    block = 'Cathode Electrolyte Anode'
  [../]
  [./bndliflux]
    order = CONSTANT
    family = MONOMIAL
    block = 'Cathode Electrolyte Anode'
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
  [./bnd_li_flux]
    type = DiffusionFluxNormalToBoundaryAux
    variable = bndliflux
    boundary = 'Electrolyte_bottom Cathode_top Electrolyte_top'
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
    block = 'Cathode Electrolyte Anode'
  [../]

  [./li_ion_flux_y]
    type = ADDiffusionFluxAux
    variable = flux_y
    component = y
    diffusion_variable = V
    diffusivity = thermal_conductivity
    block = 'Cathode Electrolyte Anode'
  [../]
  [./li_ion_flux_z]
    type = ADDiffusionFluxAux
    variable = flux_z
    component = z
    diffusion_variable = V
    diffusivity = thermal_conductivity
    block = 'Cathode Electrolyte Anode'
  [../]
[]


[Kernels]
  [./li_ion_V]
    type = ADHeatConduction
    variable = V
    use_displaced_mesh = false
    block = 'Cathode Electrolyte Anode'
  [../]
[]

[Materials]
  [./thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 100
    block = 'Cathode Anode'
  [../]
  [./thermal_conductivity2]
    type = ADHeatConductionMaterial
    thermal_conductivity = 0.1
    block = 'Electrolyte'
  [../]

  [./equilibrium_potential]
    type = ADComputeEquilibriumPotential
    R = ${gas_constant}
    faraday = ${faraday}
    temperature = ${temperature}
    cref = ${c_ref}
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
    value = ${current_density}
    variable = V
    extra_vector_tags = 'ref'
  [../]

  [./OV]
    type = ADDirichletBC
    boundary = 'Anode_top'
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
  [./bot_conc]
    type = SideAverageValue
    boundary = 'Cathode_bottom'
    variable = li_metal_conc
  [../]
  [./cathode_current]
    type = ADSideFluxIntegral
    variable = V
    diffusivity = thermal_conductivity
    boundary = 'Cathode_bottom'
  [../]
  [./anode_current]
    type = ADSideFluxIntegral
    variable = V
    diffusivity = thermal_conductivity
    boundary = 'Anode_top'
  [../]
[]

[UserObjects]
  [./vcutoff]
    type = Terminator
    expression = 'Voltage_Cathode < ${v_cutoff}'
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
  # compute_scaling_once = false
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1     NONZERO               1e-20               '
  dt = 200
  # num_steps = 20
  nl_max_its = 15
  nl_abs_tol = 1e-11
  nl_rel_tol = 1e-6
  dtmax = 200
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 15
    # timestep_limiting_postprocessor = matl_ts_min
  [../]
  end_time = 50000
  # end_time =

  # num_steps =   3
  snesmf_reuse_base = true
  scaling_group_variables = 'V V_cathode V_anode'

[]
[Problem]
  type = ReferenceResidualProblem
  # solution_variables = 'ux uy normal_lm thermal_lm li_ion_V thermal_lm2 li_metal_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'V V_cathode V_anode'
  acceptable_iterations = 2
  coord_type = RZ
[]
[Outputs]
  # exodus = true
  # csv = true
  [./csv]
    type = CSV
    file_base = csv/rct_1000_i_3.0
  [../]
  [./out]
    type = Exodus
    file_base = rst/rct_1000_i_3.0
  [../]
[]
