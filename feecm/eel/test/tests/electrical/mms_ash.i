[GlobalParams]
    energy_densities = 'E1'
[]
[Mesh]
    [Battery]
        type=GeneratedMeshGenerator
        xmax=1
        xmin=0
        nx=32
        dim=1
    []
[]
[Variables]
    [phi]
        order=FIRST
        family=LAGRANGE
    []
[]
[Kernels]
    [Laplace]
        type=RankOneDivergence
        variable = phi
        vector = i
    []
    [Source]
        type=BodyForce
        variable=phi
        function=q
    []
[]
[Functions]
    [conductivity]
        type=ParsedFunction
        expression = '1+x'
    []
    [q]
        type=ParsedFunction
        expression='2.32+exp(x)*(-7.13-3.565*x)+18.84*x'
    []
    [analyticalPhi]
        type=ParsedFunction
        expression='3.565*exp(x)-4.71*x^2+7.1*x-3.4'
    []
[]
[Materials]
    [electric_conductivity]
        type=ADGenericFunctionMaterial
        prop_names='sigma'
        prop_values='conductivity'
    []
    [electric_potential]
        type=BulkChargeTransport
        electrical_energy_density = 'E1'
        electric_potential=phi
        electric_conductivity = sigma
    []
    [Current]
        type=CurrentDensity
        current_density=i
        electric_potential=phi
    []
[]
[BCs]
    [fix]
      type = FunctionDirichletBC
      variable = phi
      boundary = 'left right'
      function = analyticalPhi
    []
  []
[Postprocessors]
    [error]
        type=ElementL2Error
        variable = phi
        function = analyticalPhi
    []
[]
[Executioner]
    type = Transient
    #solve_type = NEWTON
    solve_type = NEWTON
  
    petsc_options_iname = '-pc_type'
    petsc_options_value = 'lu'
    automatic_scaling = true
  
    num_steps = 1
[]
[Outputs]
    exodus=true
[]
