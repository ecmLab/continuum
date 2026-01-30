# ECRAM Device Simulation - Basic Template
# 3-terminal side-gated device with Li+ intercalation into 2D heterostructure
# Based on "Understanding the Mechanism of Electrochemical Intercalation into Graphene and MoS2 2D Heterostructure"

## Physical Parameters
# Gate bias
V_gate = 2.0  # V (positive for Li+ intercalation)

# PEO electrolyte properties
sigma_e = 1e-5    # S/m (ionic conductivity in 1wt%LiClO4-PEO electrolyte)
c_Li = 110        # mol/m^3 (Li+ concentration in electrolyte, LiClO4 is very dilute)
D_Li = 2.4e-14    # m^2/s (Li+ diffusivity in electrolyte, based on Nernst-Einstein relation)

# 2D material properties
sigma_2d = 1e6    # S/m (electronic conductivity of graphene)
c_max = 30000     # mol/m^3 (maximum Li+ intercalation capacity in C, LiC6)
#sigma_2d = 10    # S/m (electronic conductivity of MoS2)
#c_max = 60000     # mol/m^3 (maximum Li+ intercalation capacity in MoS2, Li1MoS2)

# Intercalation kinetics
i0 = 1e-6         # A/m^2 (exchange current density for intercalation)
alpha = 0.5       # charge transfer coefficient
U_eq = 0.2        # V (equilibrium potential for Li+ intercalation)

# Physical constants
R = 8.314         # J/(mol·K)
T = 300           # K
F = 96485         # C/mol

# Mechanical properties
E_2d = 1e12       # Pa (Young's modulus of 2D material)
nu = 0.3          # Poisson's ratio
Omega = 13e-6     # m^3/mol (molar volume of Li)

# Geometry (simplified 2D cross-section)
L_gate = 1e-6     # m (gate length)
L_channel = 10e-6 # m (channel length)
t_elyte = 100e-9  # m (electrolyte thickness)
t_2d = 1e-9       # m (2D material thickness)

[GlobalParams]
  energy_densities = 'psi_c psi_m q m'
  deformation_gradient = F
  mechanical_deformation_gradient = Fm
  eigen_deformation_gradient = Fg
  swelling_deformation_gradient = Fs
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [base]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${fparse L_gate + L_channel}
    ymin = 0
    ymax = ${t_elyte}
    nx = 50
    ny = 20
  []

  # Define gate region
  [gate_region]
    type = SubdomainBoundingBoxGenerator
    input = base
    block_id = 1
    block_name = gate
    bottom_left = '0 0 0'
    top_right = '${L_gate} ${t_elyte} 0'
  []

  # Define channel region
  [channel_region]
    type = SubdomainBoundingBoxGenerator
    input = gate_region
    block_id = 2
    block_name = channel
    bottom_left = '${L_gate} 0 0'
    top_right = '${fparse L_gate + L_channel} ${t_elyte} 0'
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
  [disp_x]   # x-displacement [m]
  []
  [disp_y]   # y-displacement [m]
  []
[]

[AuxVariables]
  [T]
    initial_condition = ${T}
  []
  [c_ref]
    initial_condition = ${c_Li}
  []
[]

[ICs]
  # Initial Li+ concentration in electrolyte
  [c_init]
    type = ConstantIC
    variable = c
    value = ${c_Li}
    block = 'gate channel'
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

  # Mechanical equilibrium
  [momentum_x]
    type = RankTwoDivergence
    variable = disp_x
    component = 0
    tensor = pk1
    factor = -1
  []
  [momentum_y]
    type = RankTwoDivergence
    variable = disp_y
    component = 1
    tensor = pk1
    factor = -1
  []
[]

[InterfaceKernels]
  # Li+ intercalation current at electrolyte/2D material interface
  [intercalation_current]
    type = MaterialInterfaceNeumannBC
    variable = Phi
    neighbor_var = Phi
    prop = i_intercalation
    factor = -1
    boundary = 'gate_channel'
  []
[]

[BCs]
  # Gate bias
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

  # Mechanical constraints
  [fix_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [fix_left]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
[]

[Materials]
  # Electrical transport
  [conductivity]
    type = ADPiecewiseConstantByBlockMaterial
    prop_name = 'sigma'
    subdomain_to_prop_value = 'gate ${sigma_e} channel ${sigma_2d}'
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

  # Chemical energy
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
    boundary = 'gate_channel'
  []

  [equilibrium_potential]
    type = ADGenericConstantMaterial
    prop_names = 'U_intercalation'
    prop_values = '${U_eq}'
    boundary = 'gate_channel'
  []

  # Mechanical properties
  [elastic_constants]
    type = ADGenericConstantMaterial
    prop_names = 'lambda G'
    prop_values = '${fparse E_2d*nu/(1+nu)/(1-2*nu)} ${fparse E_2d/2/(1+nu)}'
    block = 'channel'
  []

  [elastic_constants_electrolyte]
    type = ADGenericConstantMaterial
    prop_names = 'lambda G'
    prop_values = '1e6 1e6'  # Much softer than 2D material
    block = 'gate'
  []

  # Intercalation-induced swelling
  [swelling_coefficient_channel]
    type = ADGenericConstantMaterial
    prop_names = 'beta'
    prop_values = '0.1'  # Swelling coefficient for 2D materials
    block = 'channel'
  []
  
  [swelling_coefficient_gate]
    type = ADGenericConstantMaterial
    prop_names = 'beta'
    prop_values = '0.0'  # No swelling in electrolyte
    block = 'gate'
  []

  [swelling]
    type = SwellingDeformationGradient
    swelling_deformation_gradient = Fs
    concentration = c
    reference_concentration = c_ref
    molar_volume = ${Omega}
    swelling_coefficient = beta
  []

  [mechanical_gradient]
    type = MechanicalDeformationGradient
    mechanical_deformation_gradient = Fm
    displacements = 'disp_x disp_y'
  []
  

  [elastic_energy]
    type = NeoHookeanSolid
    elastic_energy_density = psi_m
    lambda = lambda
    shear_modulus = G
    concentration = c
    temperature = T
  []

  [stress]
    type = FirstPiolaKirchhoffStress
    first_piola_kirchhoff_stress = pk1
    deformation_gradient_rate = 'dot(F)'
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

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 15

  [TimeStepper]
    type = ConstantDT
    dt = 0.1
  []

  end_time = 10
[]

[Outputs]
  csv = true
  exodus = true
  print_linear_residuals = false
  file_base = rst/ecram_basic_full_out
[]
