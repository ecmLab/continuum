[GlobalParams]
    energy_densities='E chi'
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
    []
    [T]
    []
[]
[Kernels]
    [chargBalance]
        type=RankOneDivergence
        variable=phi
        vector=i
    []
    [SourceCharge]
        type=BodyForce
        variable=phi
        function=q
    []
    [EnergyBalance_1st_term]
        type=RankOneDivergence
        variable=T
        vector=heat_flux
    []
    [EnergyBalance_2nd_term]
        type=MaterialSource
        variable=T
        prop=r_int
        coefficient=-1
        #will comeback
    []
[]
[AuxVariables]
    [anaPhi]
        order=FIRST
        family=LAGRANGE
    []
[]
[AuxKernels]
    [anaPhiAux]
        type=FunctionAux
        variable=anaPhi
        function=analyticalPhi
    []
[]
[Functions]
    [ele_conductivity]
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
    [thermal_conductivity]
        type=ParsedFunction
        expression = '-2.6*x+3'
    []
[]
[Materials]
    [electric_conductivity]
        type=ADGenericFunctionMaterial
        prop_names='sigma'
        prop_values='ele_conductivity'
    []
    [thermal_conductivity]
        type=ADGenericFunctionMaterial
        prop_names='thermal_conductivity'
        prop_values='thermal_conductivity'
    []
    [ele_EnergyDensity]
        type=BulkChargeTransport
        electrical_energy_density = 'E'
        electric_potential=phi
        electric_conductivity = sigma
        temperature=T
    []
    [current_Force]
        type=CurrentDensity
        current_density=i
        electric_potential=phi
    []
    [xi_HeatConduction]
        type=FourierPotential
        thermal_energy_density=chi
        thermal_conductivity=thermal_conductivity
        temperature=T
    []
    [heat_flux]
        type=HeatFlux
        heat_flux=heat_flux
        temperature=T
    []
    [r_int]
        type=VariationalHeatSource
        heat_source=r_int
        temperature=T
    []
    [qconv]
        type = ADParsedMaterial
        property_name = qconv
        expression = 'htc*(T-T_inf)'
        coupled_variables = 'T'
        constant_names = 'htc T_inf'
        constant_expressions = '1.35 300'
        boundary = right
    []
[]
[BCs]
    [phiFixed]
        type=FunctionDirichletBC
        variable=phi
        boundary='left right'
        function='analyticalPhi'
    []
    [rightTemp]
        type = ADMatNeumannBC
        variable = T
        boundary = right
        value = -1
        boundary_material = qconv
    []
    # [rightTemp1]
    #     type = ConvectiveHeatFluxBC
    #     variable = T
    #     boundary = 'right'
    #     T_infinity = 300
    #     heat_transfer_coefficient = 1.35
    #     heat_transfer_coefficient_dT = 0
    # []
[]
[Postprocessors]
    [error]
      type = ElementL2Error
      variable = phi
      function = analyticalPhi
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


