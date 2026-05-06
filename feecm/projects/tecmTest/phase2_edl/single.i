# EDL model (binary electrolyte) with Poisson + Nernst–Planck for c_+ and c_-

[Mesh]
  [./base_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 50
    ny = 10
    xmax = 1.08e-8  # L_cell = 5 × xD, where xD ≈ 2.16e-9 m (Debye length)
    ymax = 2.16e-9  # Height = 1 × Debye length
    bias_x = 0.8     # Bias mesh toward boundaries
  [../]

  [./left_refine]
    type = SubdomainBoundingBoxGenerator
    input = base_mesh
    block_id = 1
    bottom_left = '0 0 0'
    top_right = '4.32e-9 2.16e-9 0'  # Refine within 2×xD from left boundary
  [../]

  [./right_refine]
    type = SubdomainBoundingBoxGenerator
    input = left_refine
    block_id = 2
    bottom_left = '6.48e-9 0 0'      # Refine within 2×xD from right boundary
    top_right = '1.08e-8 2.16e-9 0'
  [../]

  [./refine_blocks]
    type = RefineBlockGenerator
    input = right_refine
    block = '1 2'
    refinement = '2 1'  # Extra refinement in boundary regions
  [../]
[]

[Variables]
  # Electric potential
  [./phi]
    initial_condition = 0.0
    scaling = 1.0e+03
  [../]
  # Cation concentration c_+
  [./c_plus]
    initial_condition = 10.0   # mol/m^3 (bulk)
    scaling = 1.0e-01
  [../]
  # Anion concentration c_-
  [./c_minus]
    initial_condition = 10.0   # mol/m^3 (bulk)
    scaling = 1.0e-01
  [../]
[]

[AuxVariables]
  # Charge density stored as an Aux for output (also produced by material)
  [./charge_density]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxKernels]
  # Copy material charge density to Aux for output
  [./rho_from_mat]
    type = MaterialRealAux
    variable = charge_density
    property = charge_density
  [../]
[]

[Kernels]
  # Poisson: ∫ ε ∇φ·∇w dV − ∫ ρ w dV
  [./phi_diff]
    type = RankOneDivergence
    variable = phi
    vector = elec_flux
    factor = -1
  [../]
  [./phi_source]
    type = MaterialSource
    variable = phi
    prop = charge_density
    coefficient = -1.0
  [../]

  # cation c_+
  [./cplus_dt]
    type = TimeDerivative
    variable = c_plus
  [../]
  [./cplus_diff]
    type = RankOneDivergence
    variable = c_plus
    vector = j_diff_plus
    factor = -1
  [../]
  [./cplus_em]
    type = RankOneDivergence
    variable = c_plus
    vector = j_em_plus
    factor = -1
  [../]

  # anion c_-
  [./cminus_dt]
    type = TimeDerivative
    variable = c_minus
  [../]
  [./cminus_diff]
    type = RankOneDivergence
    variable = c_minus
    vector = j_diff_minus
    factor = -1
  [../]
  [./cminus_em]
    type = RankOneDivergence
    variable = c_minus
    vector = j_em_minus
    factor = -1
  [../]
[]

[BCs]
  # Stern-layer Robin BC for potential at electrode (left)
  [./phi_left]
    type = SternRobinBC
    variable = phi
    boundary = left
    epsilon0 = 8.854187817e-12
    stern_relative_permittivity = 10
    stern_thickness = 2.0e-10
    metal_potential = 1.0e-3  # 1 mV
  [../]
  # Reference potential at right (bulk)
  [./phi_right]
    type = DirichletBC
    variable = phi
    boundary = right
    value = 0.0
  [../]

  # Concentrations: blocking at left (natural), Dirichlet to bulk at right
  [./cplus_right]
    type = DirichletBC
    variable = c_plus
    boundary = right
    value = 10.0
  [../]
  [./cminus_right]
    type = DirichletBC
    variable = c_minus
    boundary = right
    value = 10.0
  [../]
[]

[Materials]
  [./electrolyte_props]
    type = GenericFunctionMaterial
    prop_names = 'D_plus D_minus permittivity'
    prop_values = '1.0e-9 2.0e-9 6.95e-10'
    block = '0 1 2'
  [../]
  [./edl_fluxes]
    type = EDLFluxes
    phi = phi
    c_plus = c_plus
    c_minus = c_minus
    permittivity = permittivity
    D_plus = D_plus
    D_minus = D_minus
    faraday_constant = 96485.33212
    gas_constant = 8.314462618
    temperature = 298.15
    z_plus = 1.0
    z_minus = -1.0
  [../]
[]

[Preconditioning]
  [./coupled]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 200
  dt = 1e-7
  solve_type = NEWTON

  # Conservative solver settings
  petsc_options_iname = '-pc_type -ksp_type'
  petsc_options_value = 'lu       preonly'

  # Solver tolerances for coupled system
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-8
  nl_max_its = 20

  l_tol = 1e-4
  l_max_its = 30

  # Steady state detection
  steady_state_detection = true
  steady_state_tolerance = 1e-7
[]

[Debug]
  show_var_residual_norms = true
[]

[Outputs]
  exodus = true
  csv = true
  file_base = rst/phase2_edl_binary
[]
