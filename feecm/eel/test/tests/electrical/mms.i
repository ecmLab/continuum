n = 32

[GlobalParams]
  energy_densities = 'E'
[]

[Mesh]
  [battery]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = 0
    xmax = 1
    nx = ${n}
  []
[]

[Variables]
  [Phi]
  []
[]

[Functions]
  [Phi]
    type = ParsedFunction
    expression = '3.565*exp(x)-4.71*x^2+7.1*x-3.4'
  []
  [sigma]
    type = ParsedFunction
    expression = '1+x'
  []
  [q]
    type = ParsedFunction
    expression = '2.32+exp(x)*(-7.13-3.565*x)+18.84*x'
  []
[]

[Kernels]
  [charge_balance]
    type = RankOneDivergence
    variable = Phi
    vector = i
  []
  [charge]
    type = BodyForce
    variable = Phi
    function = q
  []
[]

[BCs]
  [fix]
    type = FunctionDirichletBC
    variable = Phi
    boundary = 'left right'
    function = Phi
  []
[]

[Materials]
  [electric_constants]
    type = ADGenericFunctionMaterial
    prop_names = 'sigma'
    prop_values = 'sigma'
  []
  [charge_transport]
    type = BulkChargeTransport
    electrical_energy_density = E
    electric_potential = Phi
    electric_conductivity = sigma
  []
  [current]
    type = CurrentDensity
    current_density = i
    electric_potential = Phi
  []
[]

[Postprocessors]
  [error]
    type = ElementL2Error
    variable = Phi
    function = Phi
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true

  num_steps = 1
[]

[Outputs]
  exodus = true
[]
