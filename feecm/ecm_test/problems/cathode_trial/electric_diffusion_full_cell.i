# Test for full cell
# Geometry consists of
# 1) Cathode layer with given OCV curve = nmc_equilibrium_potential.csv
# 2) Solid Electrolyte
# 3) Anode at 0 V
# No mechanics included here for simplicity
# The goal is to just check to see of the whole electro-chemistry works
# All the GapDisplacementConductanceConstraint is just Linear Butler Volmer Kinetics
# All voltages are in mV => 3V = 3000 mV
[Mesh]
  [./Si]
    type = GeneratedMeshGenerator
    nx = 5
    ny = 2
    xmin = 0
    xmax = 40
    ymin = 0
    ymax = 26
    elem_type = QUAD4
    dim = 2
  [../]
  [./Si_1]
    type = RenameBlockGenerator
    input = Si
    old_block_id = '0'
    new_block_id = '1'
  [../]
  [./cathode_boundary]
    type = RenameBoundaryGenerator
    input = Si_1
    old_boundary_name = 'top bottom left right'
    new_boundary_id = '10 20 30 40'
  [../]
  [./Electrolyte]
    type = GeneratedMeshGenerator
    nx = 10
    ny = 5
    xmin = 0
    xmax = 130
    ymin = 26
    ymax = 226
    elem_type = QUAD4
    dim = 2
  [../]
  [./Electrolyte_1]
    type = RenameBlockGenerator
    input = Electrolyte
    old_block_id = '0'
    new_block_id = '2'
  [../]
  [./electrolyte_boundary]
    type = RenameBoundaryGenerator
    input = Electrolyte_1
    old_boundary_name = 'top bottom left right'
    new_boundary_id = '100 200 300 400'
  [../]
  [./anode]
    type = GeneratedMeshGenerator
    nx = 100
    ny = 5
    xmin = 0
    xmax = 11000
    ymin = 226
    ymax = 252
    elem_type = QUAD4
    dim = 2
  [../]
  [./Anode_1]
    type = RenameBlockGenerator
    input = anode
    old_block_id = '0'
    new_block_id = '3'
  [../]
  [./anode_boundary]
    type = RenameBoundaryGenerator
    input = Anode_1
    old_boundary_name = 'top bottom left right'
    new_boundary_id = '1000 2000 3000 4000'
  [../]

  [./all]
    type = CombinerGenerator
    inputs = 'cathode_boundary electrolyte_boundary anode_boundary'
  [../]
  [./primary_boundary_cathode]
    type = LowerDBlockFromSidesetGenerator
    input = all
    sidesets = 10
    new_block_name = 'cathode_s_block'
  [../]
  [./secondary_boundary_cathode]
    type = LowerDBlockFromSidesetGenerator
    input = primary_boundary_cathode
    sidesets = 200
    new_block_name = 'cathode_p_block'
  [../]
  [./secondary_anode_block]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_boundary_cathode
    sidesets = 2000
    new_block_name = 'anode_s_block'
  [../]
  [./primary_anode_block]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_anode_block
    sidesets = 100
    new_block_name = 'anode_p_block'
  [../]
  [./rename]
    type = RenameBlockGenerator
    input = primary_anode_block
    old_block_id = '1 2 3'
    new_block_name = 'Cathode Electrolyte Anode'
  [../]
  [./rename_boundary]
    type = RenameBoundaryGenerator
    input = rename
    old_boundary_id = '10 20 30 40 100 200 300 400'
    new_boundary_name = 'Cathode_top Cathode_bottom Cathode_left Cathode_right Electrolyte_top Electrolyte_bottom Electrolyte_left Electrolyte_right'
  [../]
  [./rename_boundary2]
    type = RenameBoundaryGenerator
    input = rename_boundary
    new_boundary_name = 'Anode_top Anode_bottom Anode_left Anode_right'
    old_boundary_id = '1000 2000 3000 4000'
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
    k = 1e-2
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
    k = 1e-6
    # use_displaced_mesh = true
    displacements = 'ux uy'
    compute_lm_residuals = true
    include_equilibrium_potential = false
  [../]
  # [./conc_constraint]
  #   type = ScaledBCConstraint
  #   variable = V_anode
  #   primary_variable = V
  #   secondary_variable = li_metal_conc
  #   secondary_boundary =  'Anode_bottom'
  #   secondary_subdomain = 'anode_s_block'
  #   primary_boundary = 'Electrolyte_top'
  #   primary_subdomain = 'anode_p_block'
  #
  #   primary = false
  #   scale = -1.036428e-2
  #   use_displaced_mesh = false
  # [../]
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
    block = 'Cathode Anode'
  [../]
  [./thermal_conductivity2]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1
    block = 'Electrolyte'
  [../]
  [./diffusivity_Li]
    type = ADGenericConstantMaterial
    prop_names = 'diffusivity'
    prop_values = '0.5'
    block = 'Cathode'
  [../]
  [./equilibrium_potential]
    type = ADComputeEquilibriumPotential
    R = 8.31446
    faraday = 96.4853329
    temperature = 333
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
  [./conc]
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
  nl_abs_tol = 1e-4
  nl_rel_tol = 1e-6
  dtmax = 100
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 10
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 15
    # timestep_limiting_postprocessor = matl_ts_min
  [../]
  end_time = 5000
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
  exodus = true
  csv = true
[]
