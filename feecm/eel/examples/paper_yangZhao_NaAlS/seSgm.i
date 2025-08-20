## For SE sigma comparison, only change input models but keep others the same
# sgm values are 0.01 for NPS or 0.004 for NAlS
seSgm = 0.004

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
dt = '${fparse t0/500}'      # Inital time step, to be adapted in TimeStepper, in s
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
Omega = 333  # Molar volume of Na15Sn4 from Material project, in mm^3/mmol
beta_a = 1  # Swelling coefficient of Na15Sn4, dimensionless
beta_e = 0.1  # Swelling coefficient of SE, dimensionless
beta_c = 0.00001  # Swelling coefficient of NMC, dimensionless
# Transport parameters for Na-ions
sigma_e = ${seSgm}   # Ionic conductivity of NPS, in mS/mm
sigma_a = 0.1   #??? Why would there is conductivity in the anode?? in mS/mm
sigma_c = 0.1   #??? What this mean?? assuming fast. in mS/mm
# Transport parameters for Na-atoms
D_a = 5e-7   # Na diffusivity in Na15Sn4 is 1.4e-8 mm^2/s
D_e = 1e-2   #??? Why would there is diffusivity in the electrolyte?? mm^2/s
D_c = 1e-2  # Na diffusivity in Sulfur cathode, assuming fast. mm^2/s
# Exchange current densities for the B-V reactions
#i0_a = 5e-4 # From reference paper, in mA/mm^2
i0_a = 5e-2 # From reference paper, in mA/mm^2
i0_c = 1e-1 # assuming fast. mA/mm^2
# Mechanical properties
E_a = 1.34e2 # Young's modulus of Na15Sn4, in MPa
E_e = 1.53e2 # Young's modulus of NPS, in MPa
#E_e = 1.53e3 # Young's modulus of NAlS, in MPa
E_c = 100e3  # Young's modulus of NMC, in MPa
nu_a = 0.34  # Poisson ratio of pure Na metal is 0.34
nu_e = 0.39  # From Material Project
nu_c = 0.3   # Poisson ratio of NMC

## Model geometry
#la = 0.02  # Length of Anode, in mm
#le = 0.04  # Length of SE, in mm
#lc = 0.02  # Length of cathode, in mm

## Penalty factors for solving varialbles
#c_penalty = 1
u_penalty = 1e8

[GlobalParams]
  energy_densities = 'dot(psi_m) dot(psi_c) q zeta m'
  deformation_gradient = F
  mechanical_deformation_gradient = Fm
  eigen_deformation_gradient = Fg
  swelling_deformation_gradient = Fs
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [battery]
    type = FileMeshGenerator
    file = 'data/oneDefect1um_oneSE.msh'
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
  [disp_x]
  []
  [disp_y]
  []
[]

[AuxVariables]
  [c_ref]
  []
  [T_ref]
    initial_condition = ${T0}
  []
  [pk1_xx]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = ADRankTwoAux
      rank_two_tensor = pk1
      variable = pk1_xx
      index_i = 0
      index_j = 0
    []
  []
  [stress]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = ADRankTwoScalarAux
      rank_two_tensor = pk1
      scalar_type = Hydrostatic
      execute_on = 'INITIAL TIMESTEP_END'
    []
  []
  [Js]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = ADRankTwoScalarAux
      rank_two_tensor = Fs
      scalar_type = ThirdInvariant
      execute_on = 'INITIAL TIMESTEP_END'
    []
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
  [fix_x]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = 'left right'
  []
  [fix_y]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = 'bottom'
  []
[]

[Constraints]
  [ev_y]
    type = EqualValueBoundaryConstraint
    variable = disp_y
    penalty = ${u_penalty}
    secondary = top
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
    type = ADPiecewiseConstantByBlockMaterial
    prop_name = 'beta'
    subdomain_to_prop_value = 'anode ${beta_a} elyte ${beta_e} cathode ${beta_c}'
#    prop_values = '${beta}'
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
    non_swelling_pressure = p
    output_properties = 'p'
    outputs = exodus
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
  [pressure]
    type = SideValueSampler
    variable = 'stress'
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
  end_time = ${endT}
[]

[Outputs]
  file_base = rst/oneDefect1um_oneSE_S04
  [exodus]
    type = Exodus
  []
  [csv]
    type = CSV
    execute_postprocessors_on='TIMESTEP_END'
    execute_vector_postprocessors_on='final'
  []
[]
