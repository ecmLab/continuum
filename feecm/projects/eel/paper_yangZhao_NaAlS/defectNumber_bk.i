## For defect number comparison, only change input models but keep others the same
# mode name is oneDefect_oneSE or twoDefect_oneSE
mdName = 'twoDefect1um_oneSE'

## Universal Constants
R = 8.3145 #mJ/mmol/K
T0 = 300 #K
F = 96485 #mC/mmol

## Experimental parameters
I = 0.00001   # Total cuurent, in mA
width = 0.01  # model width, in mm; assuming the model depth is 1mm
in = '${fparse I/width}'  # the current density, in mA/mm^2
C  =   1e-4   # Total capacity = in * t0, in mAh/mm^2
t0 = '${fparse abs(C/in)}'  # total charging time, in hours
dt = '${fparse t0/5000}' # Inital timestep for easy convergence, to be adapted in TimeStepper, in s
rmpT = '${fparse t0*10}'
endT = '${fparse t0*3600}'

## Material parameters
# concentrations range of pristine Materials
cmin_a = 1e-4 # minimal Na concentration in anode, mmol/mm^3
cmax_a = 4e-2 # maximal Na concentration in anode, mmol/mm^3. Pure Na is 4.2e-2
c_e    = 5e-4 #??? Na concentration in NPS, in mmol/mm^3
cmax_c = 1e-2 # maximal Na concentration in cathode, mmol/mm^3. Full Na2S is 2.4e-2
c_ref_entropy = 4.2e-2 # reference Na concentration, use Na metal
phi_eq_e = 1.12 # Equilibrium potential of SE vs Na metal, in V
# Transport parameters for Na-ions
sigma_e = 0.01   # Ionic conductivity of NPS, in mS/mm
sigma_a = 0.1   #??? Why would there is conductivity in the anode?? in mS/mm
sigma_c = 0.1   #??? What this mean?? assuming fast. in mS/mm
# Transport parameters for Na-atoms
D_a = 5e-6   # Na diffusivity in Na15Sn4 is 1.4e-8 mm^2/s
D_e = 1e-2   #??? Why would there is diffusivity in the electrolyte?? mm^2/s
D_c = 1e-2  # Na diffusivity in Sulfur cathode, assuming fast. mm^2/s
# Exchange current densities for the B-V reactions
#i0_a = 5e-4 # From reference paper, in mA/mm^2
i0_a = 5e-2 # From reference paper, in mA/mm^2
i0_c = 1e-1 # assuming fast. mA/mm^2

[GlobalParams]
  energy_densities = 'dot(psi_c) q zeta m'
[]

[Mesh]
  [battery]
    type = FileMeshGenerator
    file = 'data/${mdName}.msh'
  []
  [interfaces]
    type = BreakMeshByBlockGenerator
    input = battery
    add_interface_on_two_sides = true
    split_interface = true
  []
[]

[Variables]
  [Phi]
  []
  [c]
  []
[]

[AuxVariables]
  [c_ref]
  []
  [T_ref]
    initial_condition = ${T0}
  []
[]

[ICs]
  [c_a]
    type = ConstantIC
    variable = c
    value = ${cmin_a}
    block = 'anode'
  []
  [c_e]
    type = ConstantIC
    variable = c
    value = ${c_e}
    block = 'elyte'
  []
  [c_c]
    type = ConstantIC
    variable = c
    value = ${cmax_c}
    block = 'cathode'
  []
  [c_ref_a]
    type = ConstantIC
    variable = c_ref
    value = ${cmin_a}
    block = 'anode'
  []
  [c_ref_e]
    type = ConstantIC
    variable = c_ref
    value = ${c_e}
    block = 'elyte'
  []
  [c_ref_c]
    type = ConstantIC
    variable = c_ref
    value = ${cmax_c}
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
[]

[InterfaceKernels]
  [negative_current]
    type = MaterialInterfaceNeumannBC
    variable = Phi
    neighbor_var = Phi
    prop = ie
    factor = 1
    boundary = 'anode_elyte cathode_elyte'
  []
  [negative_mass]
    type = MaterialInterfaceNeumannBC
    variable = c
    neighbor_var = c
    prop = je
    factor = 1
    boundary = 'anode_elyte cathode_elyte'
  []
[]

[Functions]
 [in]
   type = PiecewiseLinear
   x = '0  ${rmpT}  ${endT}'  # Time points: start, ramp_end, simulation_end
   y = '0  ${in}    ${in}'  # Values: initial, ramp_target, constant_value
 []
[]

[BCs]
  [current]
    type = FunctionNeumannBC
    variable = Phi
    boundary = right
    function = in
  []
  [potential]
    type = DirichletBC
    variable = Phi
    boundary = left
    value = 0
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
  []
  [current_density]
    type = CurrentDensity
    current_density = i
    electric_potential = Phi
    output_properties = i
    outputs = exodus
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
    temperature = T_ref
    reference_concentration = ${c_ref_entropy}
    reference_chemical_potential=0.0
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
    outputs = exodus
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

  # Redox
  [ramp]
    type = ADGenericFunctionMaterial
    prop_names = 'ramp'
    prop_values = 'if(t<${t0},t/${t0},1)'
    #prop_values = 0
  []
  [OCP_anode_graphite]
    type = ADParsedMaterial
    f_name = U
#    function = 'x:=c/${cmax_a}; (2.77e-4*x^2-0.0069*x+0.0785)*ramp'
    function = 'x:=c/${cmax_a}; (2.77e-4*x^2-0.0069*x+0.0785)-${phi_eq_e}'
    args = c
    material_property_names = 'ramp'
    block = 'anode'
  []
  [OCP_cathode_NMC111]
    type = ADParsedMaterial
    f_name = U
#    function = 'x:=c/${cmax_c}; (6.0826-6.9922*x+7.1062*x^2-5.4549e-5*exp(124.23*x-114.2593)-2.5947*x^3)*ramp'
    function = 'x:=c/${cmax_c}; 6.0826-6.9922*x+7.1062*x^2-5.4549e-5*exp(124.23*x-114.2593)-2.5947*x^3-${phi_eq_e}'
    args = c
    material_property_names = 'ramp'
    block = 'cathode'
  []
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
    temperature = ${T0}
    open_circuit_potential = U
    boundary = 'anode_elyte'
    output_properties = ie
    outputs = exodus
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
    temperature = ${T0}
    open_circuit_potential = U
    boundary = 'cathode_elyte'
    output_properties = ie
    outputs = exodus
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
    boundary = left
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
[]

[VectorPostprocessors]
  [deposition_current]
    type = SideValueSampler
    variable = 'ie'
    boundary = anode_elyte
    sort_by = y
    execute_on = 'final'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true
  line_search = none

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-10
  nl_max_its = 20

  [Predictor]
    type = SimplePredictor
    scale = 1
    skip_after_failed_timestep = true
  []
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = ${dt}
    optimal_iterations = 20
    iteration_window = 5
    growth_factor = 1.5
    cutback_factor = 0.5
    cutback_factor_at_failure = 0.2
    linear_iteration_ratio = 1000
  []
  end_time = '${fparse t0*3600}'
[]

[Outputs]
  file_base = rst/${mdName}
  [exodus]
    type = Exodus
  []
  [csv]
    type = CSV
    execute_postprocessors_on='TIMESTEP_END'
    execute_vector_postprocessors_on='final'
  []
[]
