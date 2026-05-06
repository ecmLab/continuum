# Phase 1 ECRAM: Capacitive behavior (V < ECW) - TECM-no-EN Framework
# Only ion diffusion and electromigration in electrolyte, no electrochemical reactions

## Physical Parameters
V_gate = 0.01  # V (reduced voltage for faster convergence)
c_Li_e = 110      # mol/m^3 (initial concentration)
c_ClO4_e = 110    # mol/m^3
z_Li = 1          # Li+ valence
z_ClO4 = -1       # ClO4- valence
R = 8.314         # J/(mol·K) (unused)
T = 300           # K
F = 96485         # C/mol
epsilon_r = 30    # relative permittivity
epsilon_0 = 8.85e-12  # F/m
epsilon = ${fparse epsilon_r * epsilon_0}

# Transport properties - increased for faster steady state
#D_Li = 1.0e-10    # m^2/s (Li+ diffusivity)
#D_ClO4 = 1.0e-10  # m^2/s (ClO4- diffusivity)
D_Li = 2.4e-12    # m^2/s (Li+ diffusivity)
D_ClO4 = 2.0e-12  # m^2/s (ClO4- diffusivity)

# Geometry
L_total = 1e-8    # m (1 μm length)
t_device = 5e-8   # m (0.5 μm height)

[Mesh]
  [base]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${L_total}
    ymin = 0
    ymax = ${t_device}
    nx = 40
    ny = 20
  []
[]

[Variables]
  [Psi]        # Poisson potential [V]
  []
  [c_Li]       # Li+ concentration [mol/m^3]
  []
  [c_ClO4]     # ClO4- concentration [mol/m^3]
  []
[]

[AuxVariables]
  [T]
    initial_condition = ${T}
  []
  [rho_e]      # Charge density [C/m^3]
  []
  [E_x]        # Electric field x-component [V/m]
    family = MONOMIAL
    order = CONSTANT
  []
  [current_density]  # Current density [A/m^2]
    family = MONOMIAL
    order = CONSTANT
  []
  [grad_c_Li_x]      # Li+ concentration gradient x-component [mol/m^4]
    family = MONOMIAL
    order = CONSTANT
  []
  [grad_c_ClO4_x]    # ClO4- concentration gradient x-component [mol/m^4]
    family = MONOMIAL
    order = CONSTANT
  []
  [grad_Psi_x]       # Potential gradient x-component [V/m]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[ICs]
  # Start with uniform concentrations
  [c_Li_init]
    type = ConstantIC
    variable = c_Li
    value = ${c_Li_e}
  []
  [c_ClO4_init]
    type = ConstantIC
    variable = c_ClO4
    value = ${c_ClO4_e}
  []
[]

[AuxKernels]
  # Charge density calculation using TECM-no-EN kernel
  [charge_density]
    type = ChargeDensity
    variable = rho_e
    concentrations = 'c_Li c_ClO4'
    valences = '${z_Li} ${z_ClO4}'
    faraday_constant = ${F}
    fixed_charge_density = 0.0
  []

  # Electric field calculation: E_x = -dPsi/dx
  [grad_Psi_x]
    type = VariableGradientComponent
    variable = grad_Psi_x
    gradient_variable = Psi
    component = x
  []
  [E_x]
    type = ParsedAux
    variable = E_x
    coupled_variables = 'grad_Psi_x'
    expression = '-grad_Psi_x'
  []

  # Concentration gradient calculations
  [grad_c_Li_x_calc]
    type = VariableGradientComponent
    variable = grad_c_Li_x
    gradient_variable = c_Li
    component = x
  []
  [grad_c_ClO4_x_calc]
    type = VariableGradientComponent
    variable = grad_c_ClO4_x
    gradient_variable = c_ClO4
    component = x
  []

  # Correct ionic current density: j = F * (z_Li * J_Li + z_ClO4 * J_ClO4)
  # where J_i = -D_i * grad_c_i + D_i * c_i * z_i * F/(RT) * E_x
  # Note: E_x = -grad_Psi_x, so electromigration term is +D_i * c_i * z_i * F/(RT) * E_x
  [current_calc]
    type = ParsedAux
    variable = current_density
    coupled_variables = 'c_Li c_ClO4 grad_c_Li_x grad_c_ClO4_x E_x'
    expression = '${F} * (${z_Li} * (-${D_Li} * grad_c_Li_x + ${fparse D_Li * z_Li * F / (8.314 * T)} * c_Li * E_x) + ${z_ClO4} * (-${D_ClO4} * grad_c_ClO4_x + ${fparse D_ClO4 * z_ClO4 * F / (8.314 * T)} * c_ClO4 * E_x))'
  []
[]

[Kernels]
  # Poisson equation for space charge
  [poisson_equation]
    type = PoissonEquation
    variable = Psi
    permittivity = epsilon
    charge_density = rho_e
  []

  # Li+ transport: time + diffusion + electromigration
  [li_time]
    type = TimeDerivative
    variable = c_Li
  []
  [li_diffusion]
    type = ADMatDiffusion
    variable = c_Li
    diffusivity = D_Li
  []
  [li_electromigration]
    type = ADNernstPlanckConvection
    variable = c_Li
    diffusivity = D_Li
    valence = ${z_Li}
    faraday_constant = ${F}
    gas_constant = ${R}
    temperature = ${T}
    potential = Psi
  []

  # ClO4- transport: time + diffusion + electromigration
  [clo4_time]
    type = TimeDerivative
    variable = c_ClO4
  []
  [clo4_diffusion]
    type = ADMatDiffusion
    variable = c_ClO4
    diffusivity = D_ClO4
  []
  [clo4_electromigration]
    type = ADNernstPlanckConvection
    variable = c_ClO4
    diffusivity = D_ClO4
    valence = ${z_ClO4}
    faraday_constant = ${F}
    gas_constant = ${R}
    temperature = ${T}
    potential = Psi
  []
[]

[BCs]
  # Applied voltage (creates electric field)
  [psi_left]
    type = DirichletBC
    variable = Psi
    boundary = left
    value = ${V_gate}
  []
  [psi_right]
    type = DirichletBC
    variable = Psi
    boundary = right
    value = 0
  []

  # Blocking electrodes: no ion flux through boundaries
  [li_no_flux_left]
    type = NeumannBC
    variable = c_Li
    boundary = left
    value = 0
  []
  [li_no_flux_right]
    type = NeumannBC
    variable = c_Li
    boundary = right
    value = 0
  []
  [clo4_no_flux_left]
    type = NeumannBC
    variable = c_ClO4
    boundary = left
    value = 0
  []
  [clo4_no_flux_right]
    type = NeumannBC
    variable = c_ClO4
    boundary = right
    value = 0
  []
[]

[Materials]
  # Electric permittivity
  [permittivity]
    type = ADGenericConstantMaterial
    prop_names = 'epsilon'
    prop_values = '${epsilon}'
  []

  # Diffusivities
  [diffusivities]
    type = ADGenericConstantMaterial
    prop_names = 'D_Li D_ClO4'
    prop_values = '${D_Li} ${D_ClO4}'
  []

[]

[Executioner]
  type = Transient
  solve_type = NEWTON

  start_time = 0
  end_time = 1000
  dtmin = 0.0001
  dtmax = 10

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    optimal_iterations = 10
    growth_factor = 1.2
    cutback_factor = 0.5
  []

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-6
  nl_max_its = 20

  automatic_scaling = true
[]

[Postprocessors]
  # Monitor concentrations at key locations
  [li_left]
    type = PointValue
    variable = c_Li
    point = '1e-9 2.5e-8 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [li_center]
    type = PointValue
    variable = c_Li
    point = '5e-9 2.5e-8 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [li_right]
    type = PointValue
    variable = c_Li
    point = '9.9e-9 2.5e-8 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Monitor space charge
  [charge_density_max]
    type = NodalExtremeValue
    variable = rho_e
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [charge_density_min]
    type = NodalExtremeValue
    variable = rho_e
    value_type = min
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Monitor electric field
  [electric_field_max]
    type = ElementExtremeValue
    variable = E_x
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Monitor total current through device
  [total_current]
    type = SideIntegralVariablePostprocessor
    variable = current_density
    boundary = right
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Monitor current magnitude for convergence assessment
  [current_magnitude]
    type = ParsedPostprocessor
    pp_names = 'total_current'
    expression = 'abs(total_current)'
    execute_on = 'TIMESTEP_END'
  []

  # Track maximum current seen so far (to assess if current is decreasing)
  [max_current_seen]
    type = ElementExtremeValue
    variable = current_density
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Track concentration variations to quantify redistribution
  [li_concentration_spread]
    type = ParsedPostprocessor
    pp_names = 'li_left li_right'
    expression = 'abs(li_left - li_right)'
    execute_on = 'TIMESTEP_END'
  []

  [clo4_left]
    type = PointValue
    variable = c_ClO4
    point = '1e-9 2.5e-8 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [clo4_right]
    type = PointValue
    variable = c_ClO4
    point = '9.9e-9 2.5e-8 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [clo4_concentration_spread]
    type = ParsedPostprocessor
    pp_names = 'clo4_left clo4_right'
    expression = 'abs(clo4_left - clo4_right)'
    execute_on = 'TIMESTEP_END'
  []

  # Monitor relative changes from initial concentrations
  [li_relative_change_left]
    type = ParsedPostprocessor
    pp_names = 'li_left'
    expression = '(li_left - 110) / 110'
    execute_on = 'TIMESTEP_END'
  []
  [li_relative_change_right]
    type = ParsedPostprocessor
    pp_names = 'li_right'
    expression = '(li_right - 110) / 110'
    execute_on = 'TIMESTEP_END'
  []
[]

[Outputs]
  exodus = true
  csv = true
  file_base = rst/ecram_capacitive
  print_linear_residuals = false
  time_step_interval = 20
[]
