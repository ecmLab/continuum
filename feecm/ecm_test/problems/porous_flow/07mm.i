# Darcy flow with a tracer that precipitates causing mineralisation and porosity changes and permeability changes
[Mesh]
  [annular]
    type = AnnularMeshGenerator
    nr = 10
    rmin = 1.0e-6
    rmax = 10e-6
    growth_r = 1.4
    nt = 4
    dmin = 0
    dmax = 90
  []
  [make3D]
    input = annular
    type = MeshExtruderGenerator
    extrusion_vector = '0 0 12e-6'
    num_layers = 3
    bottom_sideset = 'bottom'
    top_sideset = 'top'
  []
  [shift_down]
    type = TransformGenerator
    transform = TRANSLATE
    vector_value = '0 0 -6e-6'
    input = make3D
  []
  [aquifer]
    type = SubdomainBoundingBoxGenerator
    block_id = 1
    bottom_left = '0 0 -2e-6'
    top_right = '10e-6 10e-6 2e-6'
    input = shift_down
  []
  [injection_area]
    type = ParsedGenerateSideset
    combinatorial_geometry = 'x*x+y*y<1.01e-6'
    included_subdomain_ids = 1
    new_sideset_name = 'injection_area'
    input = 'aquifer'
  []
  [rename]
    type = RenameBlockGenerator
    old_block_id = '0 1'
    new_block_name = 'caps aquifer'
    input = 'injection_area'
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [porepressure]
  []
  [tracer_concentration]
  []
[]

[PorousFlowFullySaturated]
  porepressure = porepressure
  coupling_type = Hydro
  gravity = '0 0 0'
  fp = the_simple_fluid
  mass_fraction_vars = tracer_concentration
  number_aqueous_kinetic = 1
  temperature = 283.0
  stabilization = none # Note to reader: try this with other stabilization and compare the results
[]

[AuxVariables]
  [eqm_k]
    initial_condition = 0.1
  []
  [mineral_conc]
    family = MONOMIAL
    order = CONSTANT
  []
  [initial_and_reference_conc]
    initial_condition = 0
  []
  [porosity]
    family = MONOMIAL
    order = CONSTANT
  []
  [permeability]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [mineral_conc]
    type = PorousFlowPropertyAux
    property = mineral_concentration
    mineral_species = 0
    variable = mineral_conc
  []
  [porosity]
    type = PorousFlowPropertyAux
    property = porosity
    variable = porosity
  []
  [permeability]
    type = PorousFlowPropertyAux
    property = permeability
    column = 0
    row = 0
    variable = permeability
  []
[]

[Kernels]
  [precipitation_dissolution]
    type = PorousFlowPreDis
    mineral_density = 1000.0
    stoichiometry = 1
    variable = tracer_concentration
  []
[]

[BCs]
  [constant_injection_of_tracer]
    type = PorousFlowSink
    variable = tracer_concentration
    flux_function = -5E-6
    boundary = injection_area
  []
  [constant_outer_porepressure]
    type = DirichletBC
    variable = porepressure
    value = 0
    boundary = rmax
  []
[]

[Modules]
  [FluidProperties]
    [the_simple_fluid]
      type = SimpleFluidProperties
      bulk_modulus = 2E9
      viscosity = 1.0E-3
      density0 = 1000.0
    []
  []
[]

[Materials]
  [porosity_mat]
    type = PorousFlowPorosity
    porosity_zero = 0.1
    chemical = true
    initial_mineral_concentrations = initial_and_reference_conc
    reference_chemistry = initial_and_reference_conc
  []
  [permeability_aquifer]
    type = PorousFlowPermeabilityKozenyCarman
    block = aquifer
    k0 = 1E-14
    m = 2
    n = 3
    phi0 = 0.1
    poroperm_function = kozeny_carman_phi0
  []
  [permeability_caps]
    type = PorousFlowPermeabilityKozenyCarman
    block = caps
    k0 = 1E-15
    k_anisotropy = '1 0 0  0 1 0  0 0 0.1'
    m = 2
    n = 3
    phi0 = 0.1
    poroperm_function = kozeny_carman_phi0
  []
  [precipitation_dissolution_mat]
    type = PorousFlowAqueousPreDisChemistry
    reference_temperature = 283.0
    activation_energy = 1 # irrelevant because T=Tref
    equilibrium_constants = eqm_k # equilibrium tracer concentration
    kinetic_rate_constant = 1E-8
    molar_volume = 1e-5
    num_reactions = 1
    primary_activity_coefficients = 1
    primary_concentrations = tracer_concentration
    reactions = 1
    specific_reactive_surface_area = 1
  []
  [mineral_concentration]
    type = PorousFlowAqueousPreDisMineral
  []
[]

[Preconditioning]
  active = preferred_but_might_not_be_installed
  [basic]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = ' asm      lu           NONZERO                   2'
  []
  [preferred_but_might_not_be_installed]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  []
[]

[Executioner]
  type = Transient
  automatic_scaling = true
  resid_vs_jac_scaling_param = 0.5
  solve_type = Newton
  end_time = 5
  dt = 0.1
  nl_abs_tol = 1E-15
[]

[Outputs]
  exodus = true
[]
