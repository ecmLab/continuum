gamma = 0.22
delta = 1
A = ${fparse 12*gamma/delta}
k0 = ${fparse 3*gamma*delta/2}
Refpot = 0
overPot = -0.05
phaseName = Tanh_BV_OP_${overPot}_cop_pot
L_sig = 6.25
L_eta = 0.001

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
  xmin = 0
  xmax = 10
  ymin = 0
  ymax = 10
[]

[Variables]
  [op]
    order = FIRST
    family = LAGRANGE
  []
  [pot]
    order = FIRST
    family = LAGRANGE
  []
[]

[ICs]
  [IC1]
    type = FunctionIC
    variable = op
    function = opDist
  []
  [pot]
    type = FunctionIC
    variable = pot
    function = ic_func_pot
  []
[]

[Functions]
  [opDist]
    type = ParsedFunction
    expression = '0.5*(1.0-1.0*tanh((x-x0)))'
    symbol_names = 'x0'
    symbol_values = '4'
  []
  [ic_func_pot]
    type = ParsedFunction
    expression = 'overPot*(1.0-tanh((x-x0)*2))'
    symbol_names = 'overPot x0'
    symbol_values = '${overPot} 4'
  []
[]

[Kernels]
  [op_dot]
    type = TimeDerivative
    variable = op
  []
  [bulkFree]
    type = FreeEnergyDouble
    variable = op
    A = ${A}
    scale = ${L_sig}
  []
  [lapOp]
    type = PhaseFieldLaplace
    variable = op
    k0 = ${k0}
    scale = ${L_sig}
  []
  [electrodeDr]
    type = ElectrodeDrivingForce
    variable = op
    scale = ${fparse L_eta * 8.834 * 298}
    h = h_deriv
    alpha = 0.5
    beta = 0.5
    n = 1
    F = 96485.3321
    R = 8.845
    T = 298
    conc = 1
    pot = pot
    ref_pot = ${Refpot}
  []
  # [./Cond]
  #   type = Conduction
  #   variable = pot
  #   cp=op
  #   #cv =0
  #   Q_name = eff_cond
  #   QM_name=0.
  #  [../]
    [Conduction]
      type = PhaseFieldLaplace
      variable = pot
      k0 = eff_cond
      #scale = ${L_sig}
    []
    # [coupled_pos]
    #   type = CoupledTimeDerivativePhase
    #   variable = pot
    #   phi = op
    #   n=1
    #   F=9
    #   C=1
    #   scale = -0.0074
    # []
  # [coupled_pos]
  #   type = CoupledSusceptibilityTimeDerivative
  #   variable = pot
  #   v = op
  #   f_name = 0.0074
  #   args = 'op'
  # []
[]

[AuxVariables]
  [kVal]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Materials]
  [coupled_eta_function]
    type = ADDerivativeParsedMaterial
    expression = 'op^3*(6*op^2-15*op+10)'
    coupled_variables = 'op'
    property_name = h_deriv
    derivative_order = 1
  []
  [EffectiveCond]
    type = ADDerivativeParsedMaterial
    f_name = eff_cond
    args = 'op'
    material_property_names = 'sig_s sig_l'
    function = 'sig_s*op^3*(6*op^2-15*op+10)+sig_l*(1-op^3*(6*op^2-15*op+10))'
    derivative_order = 1
  []
  [constants]
    type = ADGenericConstantMaterial
    prop_names = 'sig_s sig_l'
    prop_values = '100 1.19'
  []
[]

[BCs]
  # [left_OP]
  #   type = DirichletBC
  #   variable = op
  #   boundary = left
  #   value = 1
  # []
  # [right_OP]
  #   type = DirichletBC
  #   variable = op
  #   boundary = right
  #   value = 0
  # []
  [left_pot]
    type = DirichletBC
    variable = pot
    boundary = left
    value = ${overPot}
  []
  [right_pot]
    type = DirichletBC
    variable = pot
    boundary = right
    value = 0
  []
[]

[Postprocessors]
  [dt]
    type = TimestepSize
  []
  [aveOp]
    type = ElementAverageValue
    variable = op
  []
  [kVal]
    type = ElementAverageValue
    variable = kVal
  []
  [BulkEnergy]
    type = FBulk
    op = op
    A = ${A}
    execute_on = 'timestep_end'
  []
  [Ftot]
    type = LinearCombinationPostprocessor
    pp_names = 'BulkEnergy'
    pp_coefs = '-1'
    execute_on = 'timestep_end'
  []
  [perc_change]
    type = EnergyRatePostprocessor
    postprocessor = Ftot
    execute_on = 'timestep_end'
    dt = dt
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -snes_atol -snes_rtol -ksp_atol -ksp_rtol'
    petsc_options_value = 'hypre boomeramg 100 1e-8 1e-8 1e-8 1e-8'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  automatic_scaling = true
  start_time = 0.0
  dtmin = 1e-8
  num_steps = 50
  l_max_its = 50
  nl_max_its = 20
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    iteration_window = 2
    growth_factor = 1.5
    cutback_factor = 0.5
    dt = 1e-4
  []
[]

[Outputs]
  print_linear_residuals = false
  perf_graph_live = false
  [out]
    type = Exodus
    file_base = PhaseField_A${A}_k${k0}_${phaseName}
    elemental_as_nodal = true
  []
  [outCSV]
    type = CSV
    file_base = PhaseField_A${A}_k${k0}_${phaseName}
  []
[]
