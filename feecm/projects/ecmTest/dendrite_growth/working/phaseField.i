gamma = 0.22
delta=1
A=${fparse 12*gamma/delta}
k0=${fparse 3*gamma*delta/2}
Refpot=0
overPot=-0.01
phaseName = Tanh_BV_OP_${overPot}
#initValue = 0.6
L_sig=6.25
L_eta=0.001
[Mesh]
  # generate a 2D, 25nm x 25nm mesh
  type = GeneratedMesh
  dim = 2
  #elem_type = QUAD4
  nx = 100
  ny = 100
  # nz = 0
  xmin = 0
  xmax = 10
  ymin = 0
  ymax = 10
  # zmin = 0
  # zmax = 0
[]
# [GlobalParams]
#   # len_scale = 1.0
#    op = op
#  []
 
[Variables]
  [./op]   # Mole fraction of Cr (unitless)
    order = FIRST
    family = LAGRANGE
  [../]
[]

[ICs]
  # Use a bounding box IC at equilibrium concentrations to make sure the
  # model behaves as expected.
  [./IC1]
    type = FunctionIC
     variable = op
     function = opDist
    [../]
    # [./randomIC]
    #   type = RandomIC
    #   variable = op
    #   min=0.6
    #   max=0.8        
    # [../]
    # [./constantIC]
    #   type = ConstantIC
    #   variable = op
    #   value=${initValue}      
    # [../]
    
[]
  
[Functions]
  [./concentration]
    type = ParsedFunction
    expression = 'c0 * exp(-((x - x0)^2 / (2 * sigma_x^2))) * exp(-((y - y0)^2 / (2 * sigma_y^2)))'
    symbol_names = 'c0 x0 y0 sigma_x sigma_y'
    symbol_values = '1.0 5 0 1.0 1.0'  # Adjust values as needed
  [../]
    [./opDist]
      type = ParsedFunction
      #expression = '0.5*(1-tanh(0.5*(x-x0)))+0.5' # Making 1 to 0 order parameter makes the interface stiff
      expression = '0.5*(1.0-1.0*tanh((x-x0)))'
      symbol_names = 'x0'
      symbol_values = '8'
    [../]
[]

[Kernels]
  [./op_dot]
    type = TimeDerivative
    variable = op
  [../]
  [./bulkFree]
    variable = op
    type = FreeEnergyDouble
     A = ${A}
     scale = ${L_sig}
  [../]
  [./lapOp]
    variable = op
    type = PhaseFieldLaplace
    k0 = ${k0}
    scale = ${L_sig}
  [../]
    [./electrodeDr]
      type = ElectrodeDrivingForce
      variable = op
      scale=${fparse L_eta * 8.834 *298}
      h=h_deriv
      alpha=0.5
      beta =0.5
      n=1
      F= 96485.3321
      R= 8.845
      T= 298
      conc=1
      pot=${overPot}
      ref_pot =${Refpot}
  [../]
[]
[AuxVariables]
  [./kVal]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]
[AuxKernels]
#   [./kValAux]
#     type = ADMaterialRealAux
#     property = diffusivity
#     variable = kVal
#     execute_on = timestep_end
#   [../]
# []
[]

[Materials]
  # [diffusivity]
  #   type=AnisotropicDiffusivity
  #   # prop_names=diffusivity
  #   k0=${k0}
  #   w=0.36
  #   lambda=4
  #   gradient_variable=op
  # []
  [./coupled_eta_function]
    type = ADDerivativeParsedMaterial
    expression = 'op^3*(6*op^2-15*op+10)'
    coupled_variables = 'op'
    property_name = h_deriv
    #material_property_names = 'cs cl dh:=D[h,eta]'
    derivative_order = 1
    #outputs = exodus
  [../]

  
[]




[BCs]
  # periodic BC as is usually done on phase-field models
  # [./Periodic]
  #   [./c_bcs]
  #     auto_direction = 'x y'
  #     variable = 'op'
  #   [../]
  # [../]
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
  # [top_OP]
  #   type = NeumannBC
  #   variable = op
  #   boundary = top
  #   value = 1
  # []
  # [bottom_OP]
  #   type = NeumannBC
  #   variable = op
  #   boundary = bottom
  #   value = 0
  # []


[]
[Postprocessors]
  [./dt]
     type = TimestepSize
  [../]
  [./aveOp]
      type = ElementAverageValue
      variable =op
  [../]
  [./kVal]
    type = ElementAverageValue
    variable =kVal
  [../]
  [./BulkEnergy]
    type = FBulk
     op =op
     A=${A}
     execute_on = 'timestep_end'
  [../]
  [./Ftot]
    type = LinearCombinationPostprocessor
    pp_names = 'BulkEnergy'
    pp_coefs = ' -1'
    execute_on = 'timestep_end'
  [../]

  #   ##########################################
  #   #
  #   # NOTE: Ferret output is in attojoules
  #   #
  #   ##########################################
  # [../]
  [./perc_change]
    type = EnergyRatePostprocessor
    postprocessor = Ftot
    execute_on = 'timestep_end'
    dt = dt
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -snes_atol -snes_rtol -ksp_atol -ksp_rtol'
    petsc_options_value = 'hypre boomeramg 100 1e-8 1e-8 1e-8 1e-8'
  [../]
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
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    iteration_window = 2
    growth_factor = 1.5
    cutback_factor = 0.5
    dt = 1e-4
  [../]
[]



# [UserObjects]
#   [./kill]
#    type = Terminator
#    expression = 'perc_change <= 5.0e-8'
#   [../]
# []

[Outputs]
  print_linear_residuals = false
  perf_graph_live = false
  [./out]
    type = Exodus
    file_base = PhaseField_A${A}_k${k0}_${phaseName}
    elemental_as_nodal = true
  [../]
  [./outCSV]
    type = CSV
    file_base = PhaseField_A${A}_k${k0}_${phaseName}
  [../]
[]