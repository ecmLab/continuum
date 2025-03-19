
# Pumping a species from +z to -z via stress assisted diffusion

E = 100
nu = 0.3
lambda = '${fparse E*nu/(1+nu)/(1-2*nu)}'
mu = '${fparse E/2/(1+nu)}'

D = 100
Omega = 300
c0 = 1e-3

R = 8.3145
T = 300

[GlobalParams]
  energy_densities = 'dot(psi_m) dot(psi_c) zeta'
  deformation_gradient = F
  mechanical_deformation_gradient = Fm
  eigen_deformation_gradient = Fg
  swelling_deformation_gradient = Fs
[]

[Mesh]
  [pipe]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 6
    ny = 6
    nz = 60
    zmax = 10
  []
  [right]
    type = SideSetsFromBoundingBoxGenerator
    input = pipe
    bottom_left = '0 0 5'
    top_right = '1 1 10'
    boundaries_old = 'right'
    boundary_new = 11
    block_id = 0
  []
  [top]
    type = SideSetsFromBoundingBoxGenerator
    input = right
    bottom_left = '0 0 5'
    top_right = '1 1 10'
    boundaries_old = 'top'
    boundary_new = 12
    block_id = 0
  []
[]

[Variables]
  [c]
    initial_condition = ${c0}
  []
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
[]

[AuxVariables]
  [c0]
    initial_condition = ${c0}
  []
  [T]
    initial_condition = ${T}
  []
  [mu]
    order = CONSTANT
    family = MONOMIAL
  []
  [pressure]
    order = CONSTANT
    family = MONOMIAL
  []
[]
[AuxKernels]
  [mu_c]
    type = ADMaterialRealAux
    variable = mu
    property = mu_c
  []
  [pressureAux]
    type = ADMaterialRealAux
    variable = pressure
    property = p
  []
[]
[Kernels]
  ### Chemical
  [mass_balance_time]
    type = TimeDerivative
    variable = c
  []
  [mass_balance]
    type = RankOneDivergence
    variable = c
    vector = j
  []
  ### Mechanical
  [momentum_balance_x]
    type = RankTwoDivergence
    variable = disp_x
    tensor = P
    component = 0
    factor = -1
  []
  [momentum_balance_y]
    type = RankTwoDivergence
    variable = disp_y
    tensor = P
    component = 1
    factor = -1
  []
  [momentum_balance_z]
    type = RankTwoDivergence
    variable = disp_z
    tensor = P
    component = 2
    factor = -1
  []
[]

[BCs]
  [x_fix]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  []
  [y_fix]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  []
  [z_fix]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  []
  [push_x]
    type = FunctionNeumannBC
    variable = disp_x
    boundary = 11
    function = ramp_x
  []
  [push_y]
    type = FunctionNeumannBC
    variable = disp_y
    boundary = 12
    function = ramp_y
  []
  [push_z]
    type = FunctionNeumannBC
    variable = disp_z
    boundary = front
    function = ramp_z
  []
[]

[Functions]
  [ramp_x]
    type = PiecewiseLinear
    x = '0 0.05'
    y = '0 ${fparse -0.1*E}'
  []
  [ramp_y]
    type = PiecewiseLinear
    x = '0 0.05'
    y = '0 ${fparse -0.1*E}'
  []
  [ramp_z]
    type = PiecewiseLinear
    x = '0 0.05'
    y = '0 ${fparse -0.1*E}'
  []
[]

[Materials]
  [mobility]
    type = ADParsedMaterial
    property_name = M
    expression = '${D}*c0/${R}/${T}'
    coupled_variables = 'c0'
  []
  [chemical_energy]
    type = EntropicChemicalEnergyDensity
    chemical_energy_density = psi_c
    concentration = c
    ideal_gas_constant = ${R}
    temperature = T
    reference_concentration = 1e-4
    reference_chemical_potential = 1e3
  []
  [chemical_potential]
    type = ChemicalPotential
    chemical_potential = mu_c
    concentration = c
  []
  [diffusion]
    type = MassDiffusion
    dual_chemical_energy_density = zeta
    chemical_potential = mu_c
    mobility = M
  []
  [mass_flux]
    type = MassFlux
    mass_flux = j
    chemical_potential = mu_c
    output_properties = 'j'
    outputs = 'exodus'
  []
  [mechanical_parameters]
    type = ADGenericConstantMaterial
    prop_names = 'lambda mu beta'
    prop_values = '${lambda} ${mu} 1'
  []
  [swelling]
    type = SwellingDeformationGradient
    concentration = c
    reference_concentration = c0
    molar_volume = ${Omega}
    swelling_coefficient = beta
  []
  [def_grad]
    type = MechanicalDeformationGradient
    displacements = 'disp_x disp_y disp_z'
  []
  [neohookean]
    type = NeoHookeanSolid
    elastic_energy_density = psi_m
    lambda = lambda
    shear_modulus = mu
    concentration = c
    non_swelling_pressure = p
  []
  [pk1_stress]
    type = FirstPiolaKirchhoffStress
    first_piola_kirchhoff_stress = P
    deformation_gradient_rate = dot(F)
  []
[]
[Postprocessors]
  [conc]
    type = ElementAverageValue
    variable = c
  []
  [chemicalPot]
    type = ElementAverageValue
    variable = mu
  []
[]
[Executioner]
  type = Transient
  solve_type = NEWTON

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true

  dt = 0.001
  end_time = 0.1

  nl_rel_tol = 1e-08
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  print_linear_residuals = false
[]
