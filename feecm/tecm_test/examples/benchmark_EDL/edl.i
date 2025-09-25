# Case 1: Electric Double Layer Formation in Binary Ion System
# COMSOL benchmark validation - Gouy-Chapman-Stern theory under low voltage
# Based on COMSOL diffuse_double_layer model parameters

## Physical Parameters (from COMSOL benchmark)
T0_celsius = 25      # °C
T0 = ${fparse T0_celsius + 273.15}  # K (298.15 K)

# Ion properties
DA = 1e-9            # m²/s - Cation diffusivity
DX = 1e-9            # m²/s - Anion diffusivity (same as cation)
cA_bulk = 10         # mol/m³ - Bulk cation concentration
cX_bulk = 10         # mol/m³ - Bulk anion concentration
zA = 1               # Cation valence
zX = -1              # Anion valence

# Physical constants
R = 8.314            # J/(mol·K)
F = 96485            # C/mol
epsilon0 = 8.854e-12 # F/m - Vacuum permittivity
eps_H2O = 78.5       # Relative permittivity of water
eps_S = 10           # Relative permittivity of Stern layer
epsilon = ${fparse eps_H2O * epsilon0}  # Total permittivity

# Geometric parameters
xD = ${fparse sqrt(R * T0 * eps_H2O * epsilon0 / (2 * F * F * cA_bulk))}  # Debye length
xS = 0.2e-9          # m - Stern layer thickness (0.2 nm)
L_cell = ${fparse xD * 10}  # m - Cell length (10 × Debye length)

# Applied potential (low voltage for GCS theory validity)
phiM = 0.001         # V (1 mV)

[Mesh]
  [base]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = 0
    xmax = ${L_cell}
    nx = 100
  []
  [refine_surface]
    type = RefineBlockGenerator
    input = base
    block = 0
    refinement = 2
    enable_neighbor_refinement = false
  []
[]

[Variables]
  [phi]        # Electric potential [V]
  []
  [cA]         # Cation (A+) concentration [mol/m³]
  []
  [cX]         # Anion (X-) concentration [mol/m³]
  []
[]

[AuxVariables]
  [T]
    initial_condition = ${T0}
  []
  [rho_e]      # Charge density [C/m³]
  []
  [E_x]        # Electric field x-component [V/m]
    family = MONOMIAL
    order = CONSTANT
  []
  [deltaphi]   # Electrode-OHP potential difference [V]
    family = MONOMIAL
    order = CONSTANT
  []
  [rho_surf]   # Surface charge density [C/m²]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[ICs]
  # Start with uniform bulk concentrations (electroneutral)
  [cA_init]
    type = ConstantIC
    variable = cA
    value = ${cA_bulk}
  []
  [cX_init]
    type = ConstantIC
    variable = cX
    value = ${cX_bulk}
  []
[]

[AuxKernels]
  # Charge density calculation: rho_e = F(z_A*c_A + z_X*c_X)
  [charge_density]
    type = ChargeDensity
    variable = rho_e
    concentrations = 'cA cX'
    valences = '${zA} ${zX}'
    faraday_constant = ${F}
    fixed_charge_density = 0.0
  []

  # Electric field calculation: E_x = -dφ/dx
  [electric_field]
    type = VariableGradientComponent
    variable = E_x
    gradient_variable = phi
    component = x
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Electrode-OHP potential difference: Δφ = φ_M - φ(x=0)
  [deltaphi_calc]
    type = ParsedAux
    variable = deltaphi
    coupled_variables = 'phi'
    expression = '${phiM} - phi'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Surface charge density: σ_s = ε₀ε_S Δφ/x_S (Stern layer capacitance)
  [surface_charge_calc]
    type = ParsedAux
    variable = rho_surf
    coupled_variables = 'deltaphi'
    expression = '${fparse epsilon0 * eps_S} * deltaphi / ${xS}'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Kernels]
  # Poisson equation: -∇·(ε∇φ) = F(z_A*c_A + z_X*c_X)
  [poisson_equation]
    type = PoissonEquation
    variable = phi
    permittivity = epsilon
    charge_density = rho_e
  []

  # Cation (A+) transport: ∂c_A/∂t + ∇·J_A = 0
  # where J_A = -D_A ∇c_A - (z_A F D_A)/(RT) c_A ∇φ
  [cA_time]
    type = TimeDerivative
    variable = cA
  []
  [cA_diffusion]
    type = ADMatDiffusion
    variable = cA
    diffusivity = DA
  []
  [cA_migration]
    type = ADNernstPlanckConvection
    variable = cA
    diffusivity = DA
    valence = ${zA}
    faraday_constant = ${F}
    gas_constant = ${R}
    temperature = ${T0}
    potential = phi
  []

  # Anion (X-) transport: ∂c_X/∂t + ∇·J_X = 0
  # where J_X = -D_X ∇c_X - (z_X F D_X)/(RT) c_X ∇φ
  [cX_time]
    type = TimeDerivative
    variable = cX
  []
  [cX_diffusion]
    type = ADMatDiffusion
    variable = cX
    diffusivity = DX
  []
  [cX_migration]
    type = ADNernstPlanckConvection
    variable = cX
    diffusivity = DX
    valence = ${zX}
    faraday_constant = ${F}
    gas_constant = ${R}
    temperature = ${T0}
    potential = phi
  []
[]

[BCs]
  # For stability, let's start with simple Dirichlet boundary conditions
  # Electrode boundary (x = 0): Applied potential
  [phi_electrode]
    type = DirichletBC
    variable = phi
    boundary = left
    value = ${phiM}
  []

  # Bulk solution boundary (x = L_cell): φ = 0 (ground)
  [phi_bulk]
    type = DirichletBC
    variable = phi
    boundary = right
    value = 0
  []

  # Ion boundary conditions
  # Electrode surface: no flux (blocking electrode)
  [cA_no_flux_electrode]
    type = NeumannBC
    variable = cA
    boundary = left
    value = 0
  []
  [cX_no_flux_electrode]
    type = NeumannBC
    variable = cX
    boundary = left
    value = 0
  []

  # Bulk solution: fixed concentrations
  [cA_bulk_concentration]
    type = DirichletBC
    variable = cA
    boundary = right
    value = ${cA_bulk}
  []
  [cX_bulk_concentration]
    type = DirichletBC
    variable = cX
    boundary = right
    value = ${cX_bulk}
  []
[]

[Materials]
  # Permittivity
  [permittivity]
    type = ADGenericConstantMaterial
    prop_names = 'epsilon'
    prop_values = '${epsilon}'
  []

  # Ion diffusivities
  [diffusivities]
    type = ADGenericConstantMaterial
    prop_names = 'DA DX'
    prop_values = '${DA} ${DX}'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON

  start_time = 0
  end_time = 1
  dtmin = 1e-6
  dtmax = 10

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-6
    optimal_iterations = 6
    growth_factor = 1.1
    cutback_factor = 0.8
  []

  petsc_options_iname = '-pc_type -snes_linesearch_type -snes_rtol'
  petsc_options_value = 'lu basic 1e-4'
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-6
  nl_max_its = 20

  automatic_scaling = true
[]

[Postprocessors]
  # Monitor concentrations along the domain
  [cA_surface]
    type = PointValue
    variable = cA
    point = '${fparse xD/100} 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [cA_debye]
    type = PointValue
    variable = cA
    point = '${xD} 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [cA_bulk_check]
    type = PointValue
    variable = cA
    point = '${fparse L_cell * 0.9} 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [cX_surface]
    type = PointValue
    variable = cX
    point = '${fparse xD/100} 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [cX_debye]
    type = PointValue
    variable = cX
    point = '${xD} 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [cX_bulk_check]
    type = PointValue
    variable = cX
    point = '${fparse L_cell * 0.9} 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Monitor electric potential decay
  [phi_surface]
    type = PointValue
    variable = phi
    point = '${fparse xD/100} 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [phi_debye]
    type = PointValue
    variable = phi
    point = '${xD} 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Monitor theoretical exponential decay: φ(x_D) ≈ φ_surface * exp(-1) ≈ 0.368 * φ_surface
  [exponential_decay_check]
    type = ParsedPostprocessor
    pp_names = 'phi_surface phi_debye'
    expression = 'phi_debye / (phi_surface * exp(-1))'
    execute_on = 'TIMESTEP_END'
  []

  # Monitor charge density extrema
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

  # Monitor surface charge density
  [surface_charge]
    type = PointValue
    variable = rho_surf
    point = '0 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Calculate Debye length from current conditions (for verification)
  [debye_length_calc]
    type = ParsedPostprocessor
    pp_names = 'cA_bulk_check cX_bulk_check'
    expression = 'sqrt(${fparse R * T0 * eps_H2O * epsilon0} / (2 * ${F} * ${F} * (cA_bulk_check + cX_bulk_check) / 2))'
    execute_on = 'TIMESTEP_END'
  []
[]

[Outputs]
  exodus = true
  csv = true
  file_base = edl_results
  print_linear_residuals = false
  time_step_interval = 100
  [console]
    type = Console
    max_rows = 10
    outlier_variable_norms = false
  []
[]