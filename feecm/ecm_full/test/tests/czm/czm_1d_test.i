# Simple test for cohesize zone variable small strain formulation
# The test has 2 y elements and 1 x element
# Gap is not included in the model
# Several cases are considered
# 1) Without gap
#   a) Applied Voltage Boundary condition -> Output should be analytically computed
#      current -> Linear Butler-Volmer
#   b) Applied Voltage Boundary condition -> Output should be analytically computed
#      current -> Butler-Volmer
#   c) Applied Constant current condition -> Output should be analytically computed
#      voltage difference across the interface -> Linear Butler-Volmer
#   d) Applied Constant current condition -> Output should be analytically computed
#      voltage difference across the interface -> Butler-Volmer
# Test cases 1 a) -> 1 d) but with a concentration aux variable that changes with time
# 2) With gap
# Need to repeat in 3d as well should not matter since this is a scalar variable
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
  [subdomain_id]
    type = ParsedSubdomainMeshGenerator
    input = msh
    combinatorial_geometry = 'x > 25'
    block_id = 1
  []
  [split]
    type = BreakMeshByBlockGenerator
    input = subdomain_id
    split_interface = true
  []
[]

[GlobalParams]
  displacements = 'disp_x'
  include_gap = false
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
  [V]
    order = SECOND
  []
[]

[ICs]
  # [disp_x]
  #   type = RandomIC
  #   max = 0.0001
  #   min = -0.0001
  #   variable = disp_x
  # []
  # # [disp_y]
  # #   type = RandomIC
  # #   max = 0.0001
  # #   min = -0.0001
  # #   variable = disp_y
  # # []
  [V]
    type = ConstantIC
    # max = 10
    value = 4000.0
    variable = V
  []

[]

[AuxVariables]
  [flux_x]
    order = FIRST
    family = MONOMIAL
  []
  [flux_y]
    order = FIRST
    family = MONOMIAL
  []
  [concentration]
    block = 1
    initial_condition = 1.01e5
  []
[]

[AuxKernels]
  [flux_x]
    type = DiffusionFluxAux
    variable = flux_x
    diffusion_variable = V
    diffusivity = thermal_conductivity
    component = x
  []
  [flux_y]
    type = DiffusionFluxAux
    variable = flux_y
    diffusion_variable = V
    diffusivity = thermal_conductivity
    component = y
  []
  [concentration]
    type = FunctionAux
    function = conc
    variable = concentration
  []
[]

[Functions]
  [conc]
    type = ParsedFunction
    value = '0.74*t* 0.5e6/3600.0 + 0.22*0.5e6'
  []
  [stretch]
    type = PiecewiseLinear
    # x = '0 1'
    # y = '0 100'
    x = '0 10000 20000 30000'
    y = '0 0.5 1.0 1.5'
  []
  [eq_pot]
    type = PiecewiseLinear
    data_file = 'nmc_equilibrium_potential.csv'
    format = columns
  []
[]

[BCs]
  [fix_x]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = left
    variable = disp_x
  []

  [stretch]
    type = FunctionDirichletBC
    boundary = right
    variable = disp_x
    use_displaced_mesh = false
    function = stretch
  []

  [fix_V_top]
    type = DirichletBC
    preset = true
    value = 0.0
    variable = V
    boundary = left
  []
  # [fix_V_bottom]
  #   type = DirichletBC
  #   value = 1.0
  #   preset = true
  #   variable = V
  #   boundary = bottom
  # []
  [flux_V_bottom]
    # type = ADFunctionNeumannBC
    type = NeumannBC
    variable = V
    value = -2.4e-1
    # function = 'if (t <=1, -1e-3*t, -1e-3)'
    boundary = right
    extra_vector_tags = 'ref'
  []
  # [v_init]
  #   type = DirichletBC
  #   variable = V
  #   value = 4000.0
  #   boundary = right
  # []
[]
# [Controls]
#   [init]
#     type = TimePeriod
#     start_time = '0'
#     end_time = '0.005'
#     enable_objects = 'BCs::v_init'
#     disable_objects = 'BCs::flux_V_bottom'
#   []
# []

[Postprocessors]
  [Interface_Voltage]
    type = SideAverageValue
    variable = V
    boundary = 'Block0_Block1'
  []
  [Cell_voltage]
    type = ElementAverageValue
    variable = V
    block = 1
  []
  [Cell_concentration]
    type = ElementAverageValue
    variable = concentration
    block = 1
  []
[]

[Modules/TensorMechanics/CohesiveZoneMaster]
  [czm_ik_012]
    boundary = 'Block0_Block1'
    # generate_output = 'traction_x traction_y traction_z normal_traction tangent_traction jump_x '
    #                   'jump_y jump_z normal_jump tangent_jump'
    # base_name = 'czm_b012'
    extra_vector_tags = 'ref'
  []
[]

[InterfaceKernels]
  [czm_v]
    type = CZMDiffusiveVariableKernelSmallStrain
    boundary = 'Block0_Block1'
    neighbor_var = V
    variable = V
    include_gap = false
    component = 0
    extra_vector_tags = 'ref'
  []
[]

[Materials]
  # cohesive materials
  [czm_elastic_incremental]
    type = SalehaniIrani3DCTractionViscosity
    boundary = 'Block0_Block1'
    maximum_normal_traction = 182
    maximum_shear_traction = 364
    normal_gap_at_maximum_normal_traction = 0.01
    tangential_gap_at_maximum_shear_traction = 0.01
    normal_viscosity = 0.0549
    tangential_viscosity = 0.054
    # base_name = 'czm_b012'
  []
  # [czm_elastic_incremental]
  #   type = PureElasticTractionSeparation
  #   boundary = 'Block0_Block1'
  #   normal_stiffness = 500
  #   tangent_stiffness = 300
  #   # base_name = 'czm_b012'
  # []
  # bulk materials
  [stress]
    type = ADComputeFiniteStrainElasticStress
  []
  [elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 200e4
    poissons_ratio = 0.3
  []
  # Cohesive material for voltage
  [czm_voltage_mat]
    type = CZMButlerVolmer
    boundary = 'Block0_Block1'
    max_separation = 0.1
    k_function = "10.0"
    surfaceType = SECONDARY
    conductanceType = RCT
    computeType = "BUTLER_VOLMER"
    # gap_conductance = 100
    variable = V
    include_equil = true
    output_properties = 'flux_global interface_variable_jump'
    eq_potential_coupled_var = concentration
    reaction_rate_function = eq_pot
    cref = 0.5e6
    include_gap = false
  []
  [czm_var_jump]
    type = CZMComputeVariableJumpSmallStrain
    variable = V
    boundary = 'Block0_Block1'
  []
  [czm_global_traction]
    type = CZMComputeGlobalFluxSmallStrain
    boundary = 'Block0_Block1'
  []
  # Voltage Material
  [V_mat_1]
    type = HeatConductionMaterial
    thermal_conductivity = 0.02
    block = 0
  []
  [V_mat_2]
    type = HeatConductionMaterial
    thermal_conductivity = 1e5
    block = 1
  []

[]

[Kernels]
  [V]
    type = HeatConduction
    # block = 'blockCeramic interLayer blockMetal'
    variable = V
    use_displaced_mesh = false
  []
[]

[Modules]
  [TensorMechanics]
    [Master]
      [all]
        strain = FINITE
        add_variables = true
        use_finite_deform_jacobian = true
        use_automatic_differentiation = true
        extra_vector_tags = 'ref'
        # generate_output = 'stress_xx stress_yy stress_zz stress_xy'
      []
    []
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
  # num_steps = 10
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
[Outputs]
  exodus = true
  sync_times = '1'
  sync_only = false
  # [ex]
  #   type = Exodus
  #   output_material_properties = true
  # []
[]