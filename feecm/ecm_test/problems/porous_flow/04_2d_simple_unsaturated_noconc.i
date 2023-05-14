# Darcy flow with heat advection and conduction, and elasticity
# Dissolution of Ag in Li is the problem we are trying to solve here
# Model consists of a single block of porous material
#   1) Single phase fluid (Li/ LiAg)
#   2) Ag dissolves in the fluid
# If we want to model Ag precipitation in the fluid then perhaps we need 2 phases
[Mesh]
  [mesh]
    # type = FileMeshGenerator
    # file = test_mesh.e
    type = GeneratedMeshGenerator
    dim = 2
    nx = 60
    xmin = 1.0e-6
    xmax = 10e-6
    bias_x = 1.05
    ny = 10
    ymin = -6e-6
    ymax = 6e-6
  []
  # [annular]
  #   type = AnnularMeshGenerator
  #   nr = 10
  #   rmin = 1.0
  #   rmax = 10
  #   growth_r = 1.4
  #   nt = 4
  #   dmin = 0
  #   dmax = 90
  # []
  # [make3D]
  #   type = MeshExtruderGenerator
  #   extrusion_vector = '0 0 12'
  #   num_layers = 3
  #   bottom_sideset = 'bottom'
  #   top_sideset = 'top'
  #   input = annular
  # []
  # [shift_down]
  #   type = TransformGenerator
  #   transform = TRANSLATE
  #   vector_value = '0 0 -6'
  #   input = make3D
  # []
  [AgC]
    type = SubdomainBoundingBoxGenerator
    block_id = 1
    bottom_left = '4.99e-6 -6.1e-6 0'
    top_right = '10.1e-6 6.1e-6 0'
    input = mesh
  []
  [injection_area]
    type = ParsedGenerateSideset
    # combinatorial_geometry = 'x*x+y*y<1.01'
    combinatorial_geometry = 'x<1.0001e-6'
    included_subdomain_ids = 1
    new_sideset_name = 'injection_area'
    input = 'AgC'
  []
  [outflow_area]
    type = ParsedGenerateSideset
    combinatorial_geometry = 'x > 9.9999e-6'
    included_subdomain_ids = 1
    new_sideset_name = 'outflow_area'
    input = injection_area
  []
  [rename]
    type = RenameBlockGenerator
    old_block_id = '0 1'
    new_block_name = 'AgC LiLayer'
    input = 'outflow_area'
  []
  # [rename_boundary]
  #   type = RenameBoundaryGenerator
  #   old_boundary = '1 2 3 4'
  #   new_boundary = 'top left bottom right'
  #   input = 'rename'
  # []

[]

[Problem]
  coord_type = RZ
[]

[GlobalParams]
  displacements = 'disp_r disp_z'
  PorousFlowDictator = dictator
  biot_coefficient = 1.0
[]

[Variables]
  [porepressure]
  []
  [temperature]
    initial_condition = 293
    # scaling = 1E-5
  []

  [disp_r]
    # scaling = 1E-5
  []
  [disp_z]
    # scaling = 1E-5
  []
  # [ag_c]
  #   # initial_condition = 0.3
  # []
[]

[PorousFlowUnsaturated]
  porepressure = porepressure

  temperature = temperature
  coupling_type = ThermoHydroMechanical
  gravity = '0 0 0'
  fp = the_simple_fluid
  eigenstrain_names = thermal_contribution
  # mass_fraction_vars = ag_c
  use_displaced_mesh = false
  # number_aqueous_kinetic = 1
  multiply_by_density = true
  # stabilization = KT
  # flux_limiter_type = superbee
  relative_permeability_exponent = 3
  relative_permeability_type = Corey
  residual_saturation = 0.1
  van_genuchten_alpha = 1e-6
  van_genuchten_m = 0.6
[]

[BCs]
  [constant_injection_porepressure]
    type = DirichletBC
    variable = porepressure
    value = 0.1E6
    boundary = right
  []
  # [injected_tracer]
  #   type = DirichletBC
  #   variable = tracer_concentration
  #   value = 0.9
  #   boundary = left
  # []
  [constant_injection_flux]
    type = PorousFlowSink
    variable = porepressure
    flux_function = 'if (t <= 36000,-66.413e-9*6.8, 664.13e-9*6.8)'
    # flux_function = -664.13e-9
    boundary = left
    fluid_phase = 0
    use_relperm = true
  []
  [constant_injection_temperature]
    type = DirichletBC
    variable = temperature
    value = 293
    boundary = left
  []

  # [outflowbc]
  #   type = PorousFlowOutflowBC
  #   boundary = 'right'
  #   flux_type = fluid
  #   variable = porepressure
  #   gravity = '0 0 0'
  #   # save_in = nodal_outflow
  # []
  # [outflowbc_conc]
  #   type = PorousFlowOutflowBC
  #   boundary = 'right'
  #   flux_type = fluid
  #   variable = ag_c
  #   gravity = '0 0 0'
  #   save_in = nodal_outflow
  # []

  [roller_tmax]
    type = DirichletBC
    variable = disp_r
    value = 0
    boundary = 'left'
  []
  [roller_tmin]
    type = DirichletBC
    variable = disp_z
    value = 0
    boundary = left
  []
  [roller_top_bottom]
    type = DirichletBC
    variable = disp_z
    value = 0
    boundary = 'top bottom'
  []
  [cavity_pressure_x]
    type = Pressure
    boundary = left
    variable = disp_r
    component = 0
    factor = 1E6
    use_displaced_mesh = false
  []
  [cavity_pressure_y]
    type = Pressure
    boundary = left
    variable = disp_z
    component = 1
    factor = 1E6
    use_displaced_mesh = false
  []
[]

[Postprocessors]
  [li_mass]
    type = PorousFlowFluidMass
    PorousFlowDictator = dictator
    fluid_component = 0
  []
  # [ag_mass]
  #   type = PorousFlowFluidMass
  #   PorousFlowDictator = dictator
  #   fluid_component = 1
  # []
[]

[AuxVariables]
  [nodal_outflow]
  []
  [eqm_k]
    initial_condition = 0.1
  []
  [mineral_conc]
    family = MONOMIAL
    order = CONSTANT
  []
  [initial_and_reference_conc]
    # initial_condition = 0.2
  []
  [porosity]
    family = MONOMIAL
    order = CONSTANT
  []
  [permeability]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_xx]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_yy]
    family = MONOMIAL
    order = CONSTANT
  []
  [hydrostatic_stress]
    family = MONOMIAL
    order = CONSTANT
  []
  [vonmises_stress]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  # [mineral_conc]
  #   type = PorousFlowPropertyAux
  #   property = mineral_concentration
  #   mineral_species = 0
  #   variable = mineral_conc
  # []
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
  [stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  []
  [stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 1
    index_j = 1
  []
  [hydrostatic_stress]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    scalar_type = Hydrostatic
    variable = hydrostatic_stress
  []
  [vonmises_stress]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    scalar_type = VonMisesStress
    variable = vonmises_stress
  []

[]

# [Kernels]
#   [precipitation_dissolution]
#     type = PorousFlowPreDis
#     mineral_density = 5000.0
#     stoichiometry = 1
#     variable = ag_c
#
#   []
# []

[Modules]
  [FluidProperties]
    [the_simple_fluid]
      type = SimpleFluidProperties
      bulk_modulus = 11E9
      viscosity = 1.0e4
      density0 = 534
      thermal_expansion = 0.0000118
      cp = 4194
      cv = 4186
      porepressure_coefficient = 0
    []
  []
[]

# [ICs]
#   [ag_conc]
#     # type = RandomIC
#     type = ConstantIC
#     value = 0.2
#     # block = AgC
#     variable = ag_c
#     min = 0.01
#     max = 0.05
#   []
#   # [mineral_conc]
#   #   type = ConstantIC
#   #   value = 0.2
#   #   block = AgC
#   #   variable = mineral_conc
#   # []
# []

[Materials]
  [porosity_Agc]
    type = PorousFlowPorosityConst
    porosity = 0.5
    chemical = false
    block = AgC
    # initial_mineral_concentrations = initial_and_reference_conc
    # reference_chemistry = initial_and_reference_conc
  []

  [porosity_LiLayer]
    type = PorousFlowPorosityConst
    porosity = 0.05
    chemical = false
    block = LiLayer
    # initial_mineral_concentrations = initial_and_reference_conc
    # reference_chemistry = initial_and_reference_conc
  []

  [biot_modulus]
    type = PorousFlowConstantBiotModulus
    solid_bulk_compliance = 2E-7
    fluid_bulk_modulus = 1E7
  []
  [permeability_AgC]
    type = PorousFlowPermeabilityConst
    # block = AgC
    permeability = '1E-12 0 0   0 1E-18 0   0 0 1E-18'
  []
  # [permeability_LiLayer]
  #   type = PorousFlowPermeabilityConst
  #   block = LiLayer
  #   permeability = '1E-22 0 0   0 1E-22 0   0 0 1E-22'
  # []

  [thermal_expansion]
    type = PorousFlowConstantThermalExpansionCoefficient
    drained_coefficient = 0.00002
    fluid_coefficient = 0.00002
  []
  [rock_internal_energy]
    type = PorousFlowMatrixInternalEnergy
    density = 2500.0
    specific_heat_capacity = 1200.0
  []
  [thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '10 0 0  0 10 0  0 0 10'
    block = 'LiLayer AgC'
  []

  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 5E9
    poissons_ratio = 0.25
  []
  [strain]
    type = ComputeSmallStrain
    eigenstrain_names = thermal_contribution
  []
  [thermal_contribution]
    type = ComputeThermalExpansionEigenstrain
    temperature = temperature
    thermal_expansion_coeff = 0.00001 # this is the linear thermal expansion coefficient
    eigenstrain_name = thermal_contribution
    stress_free_temperature = 293
  []
  [stress]
    type = ComputeLinearElasticStress
  []
  # [precipitation_dissolution_mat]
  #   type = PorousFlowAqueousPreDisChemistry
  #   reference_temperature = 283.0
  #   activation_energy = 1 # irrelevant because T=Tref
  #   equilibrium_constants = eqm_k # equilibrium tracer concentration
  #   kinetic_rate_constant = 1E-3
  #   molar_volume = 10
  #   num_reactions = 1
  #   primary_activity_coefficients = 1
  #   primary_concentrations = ag_c
  #
  #   reactions = 1
  #   specific_reactive_surface_area = 1
  #   # block = AgC
  # []
  # [precipitation_dissolution_mat_li]
  #   type = PorousFlowAqueousPreDisChemistry
  #   reference_temperature = 283.0
  #   activation_energy = 1 # irrelevant because T=Tref
  #   equilibrium_constants = eqm_k # equilibrium tracer concentration
  #   kinetic_rate_constant = 1e3
  #   molar_volume = 10
  #   num_reactions = 1
  #   primary_activity_coefficients = 2
  #   primary_concentrations = ag_c
  #
  #   reactions = 1
  #   specific_reactive_surface_area = 1
  #   block = LiLayer
  # []
  # [mineral_concentration]
  #   type = PorousFlowAqueousPreDisMineral
  #   # initial_concentrations = 0.2
  #   # block = AgC
  # []
  # [mineral_concentration1]
  #   type = PorousFlowAqueousPreDisMineral
  #   initial_concentrations = 0.0
  #   block = LiLayer
  # []

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
# [Postprocessors]
#   [max_conc]
#     type = SideAverageValue
#     variable = tracer_concentration
#     boundary = 'right'
#   []
# []
# [UserObjects]
#   [term1]
#     type = Terminator
#     expression = 'max_conc > 0.95'
#   []
# []

[Executioner]
  type = Transient
  automatic_scaling = true
  # compute_scaling_once = false
  solve_type = Newton
  # end_time = 36000
  dt = 10
  nl_abs_tol = 1E-9
  nl_rel_tol = 1E-14
  resid_vs_jac_scaling_param = 0.5
  dtmax = 100
  end_time = 500
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    optimal_iterations = 15
    cutback_factor = 0.5
    growth_factor = 1.1
  []

[]

[Outputs]
  exodus = true
  [out]
    type = Checkpoint
    interval = 5
    additional_execute_on = 'FINAL'
  []
[]
