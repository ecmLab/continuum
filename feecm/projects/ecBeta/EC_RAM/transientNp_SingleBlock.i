
D_LiIon = 2.8e2
##IonicDiffusion
sig_LiIon=30
t_s=20
delta_t= 100e-1
fieldAmp = 0.2
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 400
  # ny = 14
  xmax = 200
  #ymax = 200
[]
  
[Variables]
  [./conc]
    [./InitialCondition]
      #type = ConstantIC
      type = FunctionIC
      function = guassian
      #value = 1e-4
  [../] 
  [../]
  [./pot]
  [../]

[]
[AuxVariables]
  [./T]
    initial_condition = 298
  [../]
  [./E_field]
    order = CONSTANT
    family = MONOMIAL
  [../]

[]

[AuxKernels]
  [./E_fieldAux]
    type = FunctionAux
    variable = E_field
    function = rect_pulse
  [../]
[]
  [Functions]
    # [./potDist]
    #   type = ParsedFunction
    #   expression = 'c0*(1-tanh(2*(x-x0)))'
    #   symbol_names = 'c0 x0'
    #   symbol_values = '-0.225 20'  # Adjust values as needed
    # [../]
    [guassian]
      type = ParsedFunction
      expression = 'amp*exp(-(x-mu)^2/(2*sig^2))'
      symbol_names = 'amp mu sig'
      symbol_values = '1e-3 100 20'
    []
    [./rect_pulse]
      type = ParsedFunction
      expression = 'if(t < start_time,0, if(t<=start_time+width,amplitude,0))'
      symbol_names = 'amplitude start_time width'
      symbol_values = '${fieldAmp} ${t_s} ${delta_t}'              # Amplitude of the pulse
    [../]
    [TimeVaryField]
        type = PiecewiseLinear
        x = '0 ${t_s}'
        y = '0 ${fparse 1*fieldAmp}'
    []
  []
  
  [Kernels]
    [./concTime]
      type = TimeDerivative
      variable = conc
    [../]
    [./MassFlux]
        type = ADNernstPlanckConvection
        variable = conc
        Voltage = pot
        diffusivity= ${D_LiIon}
        zIons = -1
        scale = 1
    [../]
    [./CurrentTransport]
        type = ADMigrationDiffusion
        variable = pot
        c = conc
        #T= T
        c0=1e-4
        conductivity = ${sig_LiIon}
        scale = 1
        zIons = -1
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
    #num_steps=1
    dtmin = 1e-8
    end_time = 20
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
  [BCs]
    [ElectrodeLeft]
      type=DirichletBC
      variable=pot
      boundary='left'
      value ='0'
    []
      [RightPot]
        type=DirichletBC
        variable=pot
        boundary='right'
        value =${fieldAmp}
        #function =TimeVaryField
      []
    [Flux]
        type=NeumannBC
        variable=conc
        boundary='left right'
        value ='0.0'
    []
  
  []
  [Postprocessors]
    [potSE]
      type = ElementAverageValue
      variable = pot
    []
    [appliedField]
      type = ElementAverageValue
      variable = E_field
    []
    [Conc_LiIon]
      type = ElementAverageValue
      variable = conc
    []
  []
  [Outputs]
    [./out]
      type = Exodus
      file_base = TimeVaryField_${fieldAmp}_revOFF_Gaussian
      elemental_as_nodal = true
    [../]
    [./outCSV]
      type = CSV
      file_base = TimeVaryField_${fieldAmp}_revOFF_Gaussian
      #elemental_as_nodal = true
    [../]
  []
  