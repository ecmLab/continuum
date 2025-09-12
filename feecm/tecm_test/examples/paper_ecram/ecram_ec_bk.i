# ECRAM Device Simulation - Horizontal Simplified
# 3-terminal side-gated device with Li+ intercalation into 2D heterostructure
# Left: Gate electrode, Middle: Electrolyte, Right: 2D channel material
# Only intercalation at electrolyte/channel interface

## Physical Parameters
# Gate bias (voltage-driven like working model)
V_gate = 2.0  # V (positive for Li+ intercalation)

# PEO electrolyte properties
sigma_e = 1e-5    # S/m (ionic conductivity in 1wt%LiClO4-PEO electrolyte)
c_Li_e = 110      # mol/m^3 (Li+ concentration in electrolyte, LiClO4 is very dilute)
D_Li = 2.4e-14    # m^2/s (Li+ diffusivity in electrolyte, based on Nernst-Einstein relation)

# 2D material properties
sigma_c = 1e6    # S/m (electronic conductivity of graphene)
c_Li_c = 10       # mol/m^3 (Li+ concentration in channel initially)
#sigma_c = 10    # S/m (electronic conductivity of MoS2)
#c_max = 30000     # mol/m^3 (maximum Li+ intercalation capacity in C, LiC6)
#c_max = 60000    # mol/m^3 (maximum Li+ intercalation capacity in MoS2, Li1MoS2)

# Intercalation kinetics
i0 = 1e-4         # A/m^2 (exchange current density for intercalation)
alpha = 0.5       # charge transfer coefficient
U_eq = 0.1       # V (equilibrium potential)

# Physical constants
R = 8.314         # J/(mol·K)
T = 300           # K
F = 96485         # C/mol

# Geometry (simplified 2D cross-section)
# Gate model is ignored as the gate is just the electron source
#L_gate = 2e-7     # m (gate thickness)
L_elyte = 1e-7     # m (electrolyte thickness)
L_channel = 2e-7  # m (channel thickness - 2D material)
t_device = 5e-7    # m (device height)

[GlobalParams]
  energy_densities = 'psi_c q m'
[]

[Mesh]
  [base]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${fparse L_elyte + L_channel}
    ymin = 0
    ymax = ${t_device}
    nx = 50
    ny = 92
  []

  # Define electrolyte region
  [electrolyte_region]
    type = SubdomainBoundingBoxGenerator
    input = base
    block_id = 1
    block_name = electrolyte
    bottom_left = '0 0 0'
    top_right = '${L_elyte} ${t_device} 0'
  []

  # Define channel region (2D material)
  [channel_region]
    type = SubdomainBoundingBoxGenerator
    input = electrolyte_region
    block_id = 2
    block_name = channel
    bottom_left = '${L_elyte} 0 0'
    top_right = '${fparse L_elyte + L_channel} ${t_device} 0'
  []

  # Create interface for intercalation
  [interfaces]
    type = BreakMeshByBlockGenerator
    input = channel_region
    add_interface_on_two_sides = true
    split_interface = true
  []
[]

[Variables]
  [Phi]      # Electric potential [V]
  []
  [c]        # Li+ concentration [mol/m^3]
  []
[]

[AuxVariables]
  [T]
    initial_condition = ${T}
  []
  [c_ref]
    family = MONOMIAL
    order = CONSTANT
  []
  [i_x]
    family = MONOMIAL
    order = CONSTANT
  []
  [i_y]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[ICs]
  # Initial Li+ concentration - different for electrolyte and channel
  [c_init_electrolyte]
    type = ConstantIC
    variable = c
    value = ${c_Li_e}
    block = 'electrolyte'
  []
  [c_init_channel]
    type = ConstantIC
    variable = c
    value = ${c_Li_c}
    block = 'channel'
  []
[]

[AuxKernels]
  [c_ref_aux]
    type = ADMaterialRealAux
    variable = c_ref
    property = c_ref
    execute_on = 'INITIAL'
  []
  [i_x_aux]
    type = ADMaterialRealVectorValueAux
    variable = i_x
    property = i
    component = 0
  []
  [i_y_aux]
    type = ADMaterialRealVectorValueAux
    variable = i_y
    property = i
    component = 1
  []
[]

[Kernels]
  # Charge conservation (Gauss's law)
  [charge_balance]
    type = RankOneDivergence
    variable = Phi
    vector = i
  []

  # Mass conservation for Li+
  [mass_balance_time]
    type = TimeDerivative
    variable = c
  []
  [mass_balance_flux]
    type = RankOneDivergence
    variable = c
    vector = j_total
  []
[]

[InterfaceKernels]
  # Li+ intercalation current at electrolyte/2D material interface
  [intercalation_current]
    type = MaterialInterfaceNeumannBC
    variable = Phi
    neighbor_var = Phi
    prop = i_intercalation
    factor = 1
    boundary = 'electrolyte_channel'
  []
[]

[BCs]
  # Gate bias (voltage-driven like working model)
  [gate_bias]
    type = DirichletBC
    variable = Phi
    boundary = left
    value = ${V_gate}
  []

  # Ground reference at channel end
  [channel_ground]
    type = DirichletBC
    variable = Phi
    boundary = right
    value = 0
  []
[]

[Materials]
  # Reference concentration (block-specific)
  [reference_concentration]
    type = ADPiecewiseConstantByBlockMaterial
    prop_name = 'c_ref'
    subdomain_to_prop_value = 'electrolyte ${c_Li_e} channel ${c_Li_c}'
  []

  # Electrical transport
  [conductivity]
    type = ADPiecewiseConstantByBlockMaterial
    prop_name = 'sigma'
    subdomain_to_prop_value = 'electrolyte ${sigma_e} channel ${sigma_c}'
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

  # Li+ transport in electrolyte
  [li_diffusivity]
    type = ADGenericConstantMaterial
    prop_names = 'D_Li'
    prop_values = '${D_Li}'
  []

  [li_mobility]
    type = ADParsedMaterial
    property_name = M_Li
    coupled_variables = 'c T'
    material_property_names = 'D_Li'
    expression = 'D_Li*c/${R}/T'
  []

  [chemical_energy]
    type = EntropicChemicalEnergyDensity
    chemical_energy_density = psi_c
    concentration = c
    ideal_gas_constant = ${R}
    temperature = T
    reference_concentration = c_ref
    reference_chemical_potential = 0
  []

  [chemical_potential]
    type = ChemicalPotential
    chemical_potential = mu
    concentration = c
  []

  # Mass transport (diffusion + migration)
  [diffusion]
    type = MassDiffusion
    dual_chemical_energy_density = zeta
    chemical_potential = mu
    mobility = M_Li
  []

  [migration]
    type = Migration
    electrochemical_energy_density = m
    electric_potential = Phi
    chemical_potential = mu
    electric_conductivity = sigma
    faraday_constant = ${F}
  []

  [total_mass_flux]
    type = MassFlux
    mass_flux = j_total
    energy_densities = 'zeta m'
    chemical_potential = mu
  []

  # Li+ intercalation kinetics (Butler-Volmer)
  [intercalation_kinetics]
    type = ChargeTransferReaction
    charge_transfer_current_density = i_intercalation
    charge_transfer_mass_flux = j_intercalation
    charge_transfer_heat_flux = h_intercalation
    electric_potential = Phi
    neighbor_electric_potential = Phi
    charge_transfer_coefficient = ${alpha}
    exchange_current_density = ${i0}
    faraday_constant = ${F}
    ideal_gas_constant = ${R}
    temperature = T
    open_circuit_potential = U_intercalation
    boundary = 'electrolyte_channel'
  []

  [equilibrium_potential]
    type = ADGenericConstantMaterial
    prop_names = 'U_intercalation'
    prop_values = '${U_eq}'
    boundary = 'electrolyte_channel'
  []
[]

[Postprocessors]
  # Monitor gate voltage and current
  [gate_voltage]
    type = SideAverageValue
    variable = Phi
    boundary = left
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [gate_current]
    type = ADSideIntegralMaterialProperty
    property = i
    component = 0
    boundary = left
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Monitor Li+ concentration in channel
  [li_conc_avg]
    type = ElementAverageValue
    variable = c
    block = channel
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [li_conc_max]
    type = NodalExtremeValue
    variable = c
    value_type = max
    block = channel
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON

  # Exactly like working model
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true
  line_search = none

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
  nl_max_its = 50

  # Predictor like working model
  [Predictor]
    type = SimplePredictor
    scale = 1
    skip_after_failed_timestep = true
  []

  [TimeStepper]
    type = ConstantDT
    dt = 0.1
  []

  end_time = 1.0
[]

[Outputs]
  csv = true
  exodus = true
  print_linear_residuals = false
  file_base = rst/ecram_ec
[]
