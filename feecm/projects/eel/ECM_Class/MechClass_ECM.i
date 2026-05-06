
## Current
I = 3e-4 #mA
width = 0.03 #mm
in = '${fparse I/width}'
t0 = '${fparse 1e-2/in}'
dt = '${fparse t0/10}'

cond = "FreeBC"
#Scaling the ionic conductivity by Faraday constant to avoid numerical issue
sigma_a = ${fparse 1e-8*F} #mS/mm
sigma_e = 0.1 #mS/mm
sigma_c =  ${fparse 1e-8*F} #mS/mm

# Dimension of Cell
l0 = 0
l1 = 0.04
l2 = 0.07
l3 = 0.12


##Concentration of Li Metal, scaling down by factor of 1e-1 on max value

cmin = 7.68e-3 #mmol/mm^3
cmax = 7.68e-3 #mmol/mm^3
ce   = 3.67e-5 #mmol/mm^3

## Diffusion coefficent of Li metal 
D_a = 1e-8 #mm^2/s
D_e = 1e-3 #mm^2/s
D_c = 1e-8 #mm^2/s


## Exchange current
i0_a = 1.3e-2 #mA/mm^2
i0_c = 1.3e-2 #mA/mm^2

R = 8.3145 #mJ/mmol/K
T0 = 300 #K
F = 96485 #mC/mmol


# Mechanical Properties
E_c = 12e3
E_e = 22e3
E_a = 12e3
nu_c = 0.36
nu_e = 0.37
nu_a = 0.36

# Displacement Continuity Parameter 

u_penalty = 1e8

Omega = 60  #60
## Swelling coefficient
beta_a = 10  # Anode
beta_e = 1.2e-2 # Electrolyte

[GlobalParams]
  energy_densities = 'dot(psi_m) dot(psi_c) q zeta m chi'
  deformation_gradient = F
  mechanical_deformation_gradient = Fm
  eigen_deformation_gradient = Fg
  swelling_deformation_gradient = Fs
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [battery]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = ${l0}
    xmax = ${l3}
    ymin = 0
    ymax = ${width}
    nx = 60
    ny = 15
  []
  [cathode]
    type = SubdomainBoundingBoxGenerator
    input = battery
    block_id = 1
    block_name = cathode
    bottom_left = '${l0} 0 0'
    top_right = '${l1} ${width} 0'
  []
  [elyte]
    type = SubdomainBoundingBoxGenerator
    input = cathode
    block_id = 2
    block_name = elyte
    bottom_left = '${l1} 0 0'
    top_right = '${l2} ${width} 0'
  []
  [anode]
    type = SubdomainBoundingBoxGenerator
    input = elyte
    block_id = 3
    block_name = anode
    bottom_left = '${l2} 0 0'
    top_right = '${l3} ${width} 0'
  []
  [interfaces]
    type = BreakMeshByBlockGenerator
    input = anode
    add_interface_on_two_sides = true
    split_interface = true
  []
[]

[Variables]
  [Phi]
  []
  [c]
  []
  [disp_x]
  []
  [disp_y]
  []
[]

[AuxVariables]
  [T]
    initial_condition = ${T0}
  []
  [c_ref]
  []
  [T_ref]
    initial_condition = ${T0}
  []
  [pressure]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress]
    order = CONSTANT
    family = MONOMIAL
    []
[]
[AuxKernels]
  [pressureAux]
    type = ADMaterialRealAux
    variable = pressure
    property = p
  []
  [stressKernel]
    type = ADRankTwoScalarAux
    variable = stress
    rank_two_tensor = pk1
    scalar_type = Hydrostatic
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[ICs]
  [c_min]
    type = ConstantIC
    variable = c
    value = ${cmin}
    block = 'anode'
  []
  [c_mid]
    type = ConstantIC
    variable = c
    # value = '${fparse (cmax+cmin)/2}'
    value = '${ce}'
    block = 'elyte'
  []
  [c_max]
    type = ConstantIC
    variable = c
    value = ${cmax}
    block = 'cathode'
  []
  [c_ref_min]
    type = ConstantIC
    variable = c_ref
    value = ${cmin}
    block = 'anode'
  []
  [c_ref_mid]
    type = ConstantIC
    variable = c_ref
    #value = '${fparse (cmax+cmin)/2}'
    value = '${ce}'
    block = 'elyte'
  []
  [c_ref_max]
    type = ConstantIC
    variable = c_ref
    value = ${cmax}
    block = 'cathode'
  []
[]

[Kernels]
  # Charge balance
  [charge_balance]
    type = RankOneDivergence
    variable = Phi
    vector = i
  []
  # Mass balance
  [mass_balance_1]
    type = TimeDerivative
    variable = c
  []
  [mass_balance_2]
    type = RankOneDivergence
    variable = c
    vector = j
  []
  # Momentum balance
  [momentum_balance_x]
    type = RankTwoDivergence
    variable = disp_x
    component = 0
    tensor = pk1
    factor = -1
  []
  [momentum_balance_y]
    type = RankTwoDivergence
    variable = disp_y
    component = 1
    tensor = pk1
    factor = -1
  []
[]

[InterfaceKernels]
  [negative_current]
    type = MaterialInterfaceNeumannBC
    variable = Phi
    neighbor_var = Phi
    prop = ie
    factor = -1
    boundary = 'elyte_anode cathode_elyte'
  []
  [continuity_disp_x]
    type = InterfaceContinuity
    variable = disp_x
    neighbor_var = disp_x
    penalty = ${u_penalty}
    boundary = 'anode_elyte elyte_cathode'
  []
  [continuity_disp_y]
    type = InterfaceContinuity
    variable = disp_y
    neighbor_var = disp_y
    penalty = ${u_penalty}
    boundary = 'anode_elyte elyte_cathode'
  []
[]
[Constraints]
  [y]
    type = EqualValueBoundaryConstraint
    variable = disp_y
    penalty = ${u_penalty}
    secondary = top
  []
[]

[Functions]
  [ramp_current]
    type = PiecewiseLinear
    x = '0 ${fparse 1*t0}'
    y = '0 ${fparse 1*in}'
  []
  [pulse]
    type = PiecewiseConstant
    x = '0 300 600'
    y = '${fparse 1*in} ${fparse -1*in} ${fparse 1*in}'
  []
[]

[BCs]
  [left]
    type = FunctionNeumannBC
    variable = Phi
    boundary = left
    function = pulse 
  []
  [right]
    type = DirichletBC
    variable = Phi
    boundary = right
    value = 0
  []
  [open]
    type = OpenBC
    variable = c
    flux = jm
    boundary = 'left right'
  []
  [fix_x]
    type = DirichletBC
    variable = 'disp_x'
    value = 0
    boundary = 'left'
  []
  # [pressure_x]
  #   type = Pressure
  #   variable = 'disp_x'
  #   factor = ${pressure}
  #   boundary = 'right'
  # []
  [fix_y]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = 'bottom top'
  []
[]

[Materials]
  # Electrodynamics
  [conductivity]
    type = ADPiecewiseConstantByBlockMaterial
    prop_name = 'sigma'
    subdomain_to_prop_value = 'anode ${sigma_a} elyte ${sigma_e} cathode ${sigma_c}'
  []
  [charge_transport]
    type = BulkChargeTransport
    electrical_energy_density = q
    electric_potential = Phi
    electric_conductivity = sigma
    temperature = T
  []
  [current_density]
    type = CurrentDensity
    current_density = i
    electric_potential = Phi
  []

  # Chemical reactions
  [diffusivity]
    type = ADPiecewiseConstantByBlockMaterial
    prop_name = 'D'
    subdomain_to_prop_value = 'anode ${D_a} elyte ${D_e} cathode ${D_c}'
  []
  [mobility]
    type = ADParsedMaterial
    f_name = M
    args = 'c_ref T_ref'
    material_property_names = 'D'
    function = 'D*c_ref/${R}/T_ref'
  []
  [chemical_energy]
    type = EntropicChemicalEnergyDensity
    chemical_energy_density = psi_c
    concentration = c
    ideal_gas_constant = ${R}
    temperature = T
    reference_concentration = c_ref
    reference_chemical_potential=0
  []
  [chemical_potential]
    type = ChemicalPotential
    chemical_potential = mu
    concentration = c
  []
  [diffusion]
    type = MassDiffusion
    dual_chemical_energy_density = zeta
    chemical_potential = mu
    mobility = M
  []
  [mass_flux]
    type = MassFlux
    mass_flux = j
    chemical_potential = mu
  []

  # Migration
  [migration]
    type = Migration
    electrochemical_energy_density = m
    electric_potential = Phi
    chemical_potential = mu
    electric_conductivity = sigma
    faraday_constant = ${F}
  []
  [migration_flux]
    type = MassFlux
    mass_flux = jm
    energy_densities = 'm'
    chemical_potential = mu
  []
  [ramp]
    type = ADGenericFunctionMaterial
    prop_names = 'ramp'
    prop_values = 'if(t<${t0},t/${t0},1)'
    #prop_values = 0
  []
  [./OCP_anode_Li]
    type = ADGenericConstantMaterial
    prop_names = 'U'
    prop_values = '0'
    block = 'anode'
  [../]
  [./OCP_cathode]
    type = ADGenericConstantMaterial
    prop_names = 'U'
    prop_values = '0'
    block = 'cathode'
  [../]
  [charge_transfer_elyte_anode]
    type = ChargeTransferReaction
    charge_transfer_current_density = ie
    charge_transfer_mass_flux = je
    charge_transfer_heat_flux = he
    electric_potential = Phi
    neighbor_electric_potential = Phi
    charge_transfer_coefficient = 0.5
    exchange_current_density = ${i0_a}
    faraday_constant = ${F}
    ideal_gas_constant = ${R}
    temperature = T
    open_circuit_potential = U
    boundary = 'elyte_anode'
  []
  [charge_transfer_cathode_elyte]
    type = ChargeTransferReaction
    charge_transfer_current_density = ie
    charge_transfer_mass_flux = je
    charge_transfer_heat_flux = he
    electric_potential = Phi
    neighbor_electric_potential = Phi
    charge_transfer_coefficient = 0.5
    exchange_current_density = ${i0_c}
    faraday_constant = ${F}
    ideal_gas_constant = ${R}
    temperature = T
    open_circuit_potential = U
    boundary = 'cathode_elyte'
  []
  # Mechanical
  [swelling_coef]
    type = ADPiecewiseConstantByBlockMaterial
    prop_name = 'beta'
    subdomain_to_prop_value = 'anode ${beta_a} elyte ${beta_e} cathode ${beta_a}'
  []
  [stiffness_c]
    type = ADGenericConstantMaterial
    prop_names = 'lambda G'
    prop_values = '${fparse E_c*nu_c/(1+nu_c)/(1-2*nu_c)} ${fparse E_c/2/(1+nu_c)}'
    block = cathode
  []
  [stiffness_e]
    type = ADGenericConstantMaterial
    prop_names = 'lambda G'
    prop_values = '${fparse E_e*nu_e/(1+nu_e)/(1-2*nu_e)} ${fparse E_e/2/(1+nu_e)}'
    block = elyte
  []
  [stiffness_a]
    type = ADGenericConstantMaterial
    prop_names = 'lambda G'
    prop_values = '${fparse E_a*nu_a/(1+nu_a)/(1-2*nu_a)} ${fparse E_a/2/(1+nu_a)}'
    block = anode
  []
  [swelling]
    type = SwellingDeformationGradient
    concentration = c
    reference_concentration = c_ref
    molar_volume = ${Omega}
    swelling_coefficient = beta
  []
  [defgrad]
    type = MechanicalDeformationGradient
    displacements = 'disp_x disp_y'
  []
  [neohookean]
    type = NeoHookeanSolid
    elastic_energy_density = psi_m
    lambda = lambda
    shear_modulus = G
    concentration = c
    temperature = T_ref
    non_swelling_pressure = p
  []
  [pk1]
    type = FirstPiolaKirchhoffStress
    first_piola_kirchhoff_stress = pk1
    deformation_gradient_rate = dot(F)
  []
[]

[Postprocessors]
  [V_l]
    type = SideAverageValue
    variable = Phi
    boundary = left
    outputs = none
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [V_r]
    type = SideAverageValue
    variable = Phi
    boundary = right
    outputs = none
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [V]
    type = ParsedPostprocessor
    function = 'V_l - V_r'
    pp_names = 'V_l V_r'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [I]
    type = ADSideIntegralMaterialProperty
    property = i
    component = 0
    boundary = right
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [dt]
    type = TimestepSize
    outputs = none
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [dC]
    type = ParsedPostprocessor
    function = '-dt*I'
    pp_names = 'dt I'
    outputs = none
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [C]
    type = CumulativeValuePostprocessor
    postprocessor = dC
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_a_max]
    type = NodalExtremeValue
    variable = c
    value_type = max
    block = anode
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_c_min]
    type = NodalExtremeValue
    variable = c
    value_type = min
    block = cathode
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_a_min]
    type = NodalExtremeValue
    variable = c
    value_type = min
    block = anode
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_c_max]
    type = NodalExtremeValue
    variable = c
    value_type = max
    block = cathode
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [mass_a]
    type = ElementIntegralVariablePostprocessor
    variable = c
    block = anode
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [mass_e]
    type = ElementIntegralVariablePostprocessor
    variable = c
    block = elyte
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [mass_c]
    type = ElementIntegralVariablePostprocessor
    variable = c
    block = cathode
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [pulse_out]
    type = FunctionValuePostprocessor
    function = pulse
    execute_on = 'initial timestep_end'
  []
[]

[UserObjects]
  [kill_voltage]
   type = Terminator
   expression = 'V >= 5'
   message = 'Concentration in cathode is below the minimum allowable value.'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true
  line_search = none

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
  nl_max_its = 20

  [Predictor]
    type = SimplePredictor
    scale = 1
    skip_after_failed_timestep = true
  []
  [TimeStepper]
    type = ConstantDT
    dt = ${fparse 10*dt}
  []
  end_time = 900
[]

[Outputs]
  file_base = 'ECM_class/ECM_I_${I}__LiMetal_${cond}'
  csv = true
  exodus = true
  print_linear_residuals = false
  checkpoint = true
[]
