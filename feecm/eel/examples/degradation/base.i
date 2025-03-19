I = 5e-4 #mA
width = 0.03 #mm
in = '${fparse -I/width}'
t0 = '${fparse -1e-2/in}'
dt = '${fparse t0/50}'

sigma_a = 1e0 #mS/mm
sigma_e = 1e-1 #mS/mm
sigma_c = 1e-2 #mS/mm

l0 = 0
l1 = 0.04
l2 = 0.07
l3 = 0.12

cmin = 1e-4 #mmol/mm^3
cmax = 1e-3 #mmol/mm^3
D_a = 1e-3 #mm^2/s
D_e = 1e-4 #mm^2/s
D_c = 5e-5 #mm^2/s
mu0 = 1e3

h0 = 1e-9
hc = 1e-3
Omega_sei = 100
A_sei = 0.8
Q_sei = 5e3
rho0_sei = 5e2

R = 8.3145 #mJ/mmol/K
T0 = 300 #K
F = 96485 #mC/mmol

i0_a = 1e-4 #mA/mm^2
i0_c = 1e-4 #mA/mm^2

E_c = 1e5
E_e = 1e4
E_a = 2e5
nu_c = 0.3
nu_e = 0.25
nu_a = 0.3

Ei = 1e7
Gi = 5e6
Acr = 1e-5
Qcr = 4e3
ncr = 5
Tn0 = 10
Gc = 2.06e-5

u_penalty = 1e8

Omega = 60
beta = 1e-3
CTE = 1e-5

[GlobalParams]
  energy_densities = 'dot(psi_m) dot(psi_c) q zeta'
  displacements = 'disp_x disp_y'
  deformation_gradient = F
  mechanical_deformation_gradient = Fm
  eigen_deformation_gradient = Fg
  swelling_deformation_gradient = Fs
  thermal_deformation_gradient = Ft
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
    ny = 1
  []
  [anode]
    type = SubdomainBoundingBoxGenerator
    input = battery
    block_id = 1
    block_name = anode
    bottom_left = '${l0} 0 0'
    top_right = '${l1} ${width} 0'
  []
  [elyte]
    type = SubdomainBoundingBoxGenerator
    input = anode
    block_id = 2
    block_name = elyte
    bottom_left = '${l1} 0 0'
    top_right = '${l2} ${width} 0'
  []
  [cathode]
    type = SubdomainBoundingBoxGenerator
    input = elyte
    block_id = 3
    block_name = cathode
    bottom_left = '${l2} 0 0'
    top_right = '${l3} ${width} 0'
  []
  [interfaces]
    type = BreakMeshByBlockGenerator
    input = cathode
    add_interface_on_two_sides = true
    split_interface = true
  []
  use_displaced_mesh = false
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
  [c_ref]
  []
  [Phi0]
  []
  [T]
    initial_condition = ${T0}
  []
  [T_ref]
    initial_condition = ${T0}
  []
[]

[ICs]
  [c_ref_min]
    type = ConstantIC
    variable = c_ref
    value = ${cmin}
    block = 'anode'
  []
  [c_ref_mid]
    type = ConstantIC
    variable = c_ref
    value = '${fparse (cmax+cmin)/2}'
    block = 'elyte'
  []
  [c_ref_max]
    type = ConstantIC
    variable = c_ref
    value = ${cmax}
    block = 'cathode'
  []
[]

[AuxKernels]
  [Phi0]
    type = ParsedAux
    variable = Phi0
    expression = 'Phi'
    coupled_variables = 'Phi'
    execute_on = 'INITIAL'
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

[Modules]
  [TensorMechanics]
    [CohesiveZoneMaster]
      [interface]
        boundary = 'anode_elyte'
        strain = SMALL
        use_automatic_differentiation = true
      []
    []
  []
[]

[InterfaceKernels]
  [negative_current]
    type = MaterialInterfaceNeumannBC
    variable = Phi
    neighbor_var = Phi
    prop = ie
    factor = -1
    boundary = 'cathode_elyte'
  []
  [positive_current]
    type = MaterialInterfaceNeumannBC
    variable = Phi
    neighbor_var = Phi
    prop = ie
    boundary = 'anode_elyte'
  []
  [negative_mass]
    type = MaterialInterfaceNeumannBC
    variable = c
    neighbor_var = c
    prop = je
    factor = -1
    boundary = 'cathode_elyte'
  []
  [positive_mass]
    type = MaterialInterfaceNeumannBC
    variable = c
    neighbor_var = c
    prop = je
    boundary = 'anode_elyte'
  []
  [continuity_disp_x]
    type = InterfaceContinuity
    variable = disp_x
    neighbor_var = disp_x
    penalty = ${u_penalty}
    boundary = 'cathode_elyte'
  []
  [continuity_disp_y]
    type = InterfaceContinuity
    variable = disp_y
    neighbor_var = disp_y
    penalty = ${u_penalty}
    boundary = 'cathode_elyte'
  []
[]

[BCs]
  [pull_x]
    type = NeumannBC
    variable = disp_x
    value = 1
    boundary = 'right'
  []
  [fix_x]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = 'left'
  []
  [fix_y]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = 'bottom'
  []
[]

[Functions]
  [in]
    type = PiecewiseLinear
    x = '0 ${t0}'
    y = '0 ${in}'
  []
  [in_discharging]
    type = PiecewiseLinear
    x = '0 ${t0}'
    y = '0 ${fparse -in}'
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
    property_name = M
    coupled_variables = 'c_ref T_ref'
    material_property_names = 'D'
    expression = 'D*c_ref/${R}/T_ref'
  []
  [chemical_energy]
    type = EntropicChemicalEnergyDensity
    chemical_energy_density = psi_c
    concentration = c
    ideal_gas_constant = ${R}
    temperature = T
    reference_concentration = c_ref
    reference_chemical_potential = ${mu0}
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

  # Redox
  [OCP_anode_graphite]
    type = ADParsedMaterial
    f_name = U
    function = 'x:=c/${cmax}; (2.77e-4*x^2-0.0069*x+0.0785)'
    args = c
    material_property_names = 'ramp'
    boundary = 'anode_elyte'
  []
  [OCP_cathode_NMC111]
    type = ADParsedMaterial
    f_name = U
    function = 'x:=c/${cmax}; (6.0826-6.9922*x+7.1062*x^2-5.4549e-5*exp(124.23*x-114.2593)-2.5947*x^3)*ramp'
    args = c
    material_property_names = 'ramp'
    boundary = 'cathode_elyte'
  []
  [charge_transfer_anode_elyte]
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
    interface_resistance = rho_sei
    degradation = g
    boundary = 'anode_elyte'
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
  [SEI]
    type = SEIGrowth
    thickness = h_sei
    initial_thickness = ${h0}
    charge_transfer_mass_flux = je
    characteristic_thickness = ${hc}
    molar_volume = ${Omega_sei}
    correction = ${A_sei}
    activation_energy = ${Q_sei}
    ideal_gas_constant = ${R}
    temperature = T
    boundary = 'anode_elyte'
  []
  [SEI_resistance]
    type = ADParsedMaterial
    property_name = rho_sei
    expression = '${rho0_sei} * h_sei'
    material_property_names = 'h_sei'
    boundary = 'anode_elyte'
  []

  # Mechanical
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
  [swelling_coefficient]
    type = ADGenericConstantMaterial
    prop_names = 'beta'
    prop_values = '${beta}'
  []
  [swelling]
    type = SwellingDeformationGradient
    concentration = c
    reference_concentration = c_ref
    molar_volume = ${Omega}
    swelling_coefficient = beta
  []
  [thermal_expansion]
    type = ThermalDeformationGradient
    temperature = T
    reference_temperature = T_ref
    CTE = ${CTE}
  []
  [defgrad]
    type = MechanicalDeformationGradient
  []
  [neohookean]
    type = NeoHookeanSolid
    elastic_energy_density = psi_m
    lambda = lambda
    shear_modulus = G
    concentration = c
    temperature = T
  []
  [pk1]
    type = FirstPiolaKirchhoffStress
    first_piola_kirchhoff_stress = pk1
    deformation_gradient_rate = dot(F)
  []
  [traction]
    type = InterfaceTractionWithCreepDegradation
    normal_stiffness = ${Ei}
    tangential_stiffness = ${Gi}
    normal_traction = Tn
    activation_energy = ${Qcr}
    ideal_gas_constant = ${R}
    temperature = T
    damage = Di
    degradation = g
    creep_displacement_jump = juc
    fracture_driving_energy = psi_f
    energy_release_rate = ${Gc}
    reference_normal_traction = ${Tn0}
    creep_coefficient = ${Acr}
    creep_exponent = ${ncr}
    residual_stiffness = 1e-4
    boundary = 'anode_elyte'
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
    function = 'V_r - V_l'
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
    outputs = none
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_c_min]
    type = NodalExtremeValue
    variable = c
    value_type = min
    block = cathode
    outputs = none
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_a_min]
    type = NodalExtremeValue
    variable = c
    value_type = min
    block = anode
    outputs = none
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c_c_max]
    type = NodalExtremeValue
    variable = c
    value_type = max
    block = cathode
    outputs = none
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
  [SEI_thickness]
    type = ADSideAverageMaterialProperty
    property = h_sei
    boundary = 'anode_elyte'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Di]
    type = ADSideAverageMaterialProperty
    property = Di
    boundary = 'anode_elyte'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [psi_f]
    type = ADSideAverageMaterialProperty
    property = psi_f
    boundary = 'anode_elyte'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [dV]
    type = ChangeOverTimePostprocessor
    postprocessor = V
    outputs = none
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [target_V_dt]
    type = ParsedPostprocessor
    function = 'if(abs(dV)>0.01,${dt},1e6)'
    pp_names = 'dV'
    outputs = none
    execute_on = 'TIMESTEP_END'
  []
[]
