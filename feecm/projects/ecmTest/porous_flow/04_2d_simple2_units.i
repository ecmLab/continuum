# Darcy flow with heat advection and conduction, and elasticity
# Units ---
# Length - um
# viscosity - 1 cP -> 1e-3 Pa-s -> 1e-9 Pa-
[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    xmin = 1.0e-6
    xmax = 10e-6
    bias_x = 1.2
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
  [aquifer]
    type = SubdomainBoundingBoxGenerator
    block_id = 1
    bottom_left = '0 -2e-6 0'
    top_right = '10e-6 2e-6 0'
    input = mesh
  []
  [injection_area]
    type = ParsedGenerateSideset
    # combinatorial_geometry = 'x*x+y*y<1.01'
    combinatorial_geometry = 'x<1.0001e-6'
    included_subdomain_ids = 1
    new_sideset_name = 'injection_area'
    input = 'aquifer'
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
    new_block_name = 'caps aquifer'
    input = 'outflow_area'
  []

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
  # [tracer_concentration]
  # []
  [disp_r]
    # scaling = 1E-5
  []
  [disp_z]
    # scaling = 1E-5
  []
[]

[PorousFlowFullySaturated]
  porepressure = porepressure
  temperature = temperature
  coupling_type = ThermoHydroMechanical
  gravity = '0 0 0'
  fp = the_simple_fluid
  eigenstrain_names = thermal_contribution
  # mass_fraction_vars = tracer_concentration
  use_displaced_mesh = false
  # relative_permeability_exponent = 3
  # relative_permeability_type = Corey
  # residual_saturation = 0.9
  # van_genuchten_alpha = 1E-6
  # van_genuchten_m = 0.6
[]

[BCs]
  [constant_injection_porepressure]
    type = DirichletBC
    variable = porepressure
    value = 1E6
    boundary = right
  []
  # [injected_tracer]
  #   type = DirichletBC
  #   variable = tracer_concentration
  #   value = 0.5
  #   boundary = left
  # []
  [constant_injection_flux]
    type = PorousFlowSink
    variable = porepressure
    # flux_function = 'if (t <= 1e6,-1.0e-2, 1.0e-2)'
    flux_function = -66e-9 # kg/m^2s
    boundary = left
    fluid_phase = 0
  []
  [constant_injection_temperature]
    type = DirichletBC
    variable = temperature
    value = 313
    boundary = left
  []

  [outflowbc]
    type = PorousFlowOutflowBC
    boundary = 'right'
    flux_type = fluid
    variable = porepressure
    gravity = '0 0 0'
    # save_in = nodal_outflow
  []
  # [outflowbc_conc]
  #   type = PorousFlowOutflowBC
  #   boundary = 'right'
  #   flux_type = fluid
  #   variable = tracer_concentration
  #   gravity = '0 0 0'
  #   # save_in = nodal_outflow
  # []

  [roller_tmax]
    type = DirichletBC
    variable = disp_r
    value = 0
    boundary = 'left right'
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

[AuxVariables]
  [stress_rr]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_pp]
    family = MONOMIAL
    order = CONSTANT
  []
[]

# [AuxKernels]
#   [stress_rr]
#     type = RankTwoScalarAux
#     rank_two_tensor = stress
#     variable = stress_rr
#     scalar_type = RadialStress
#     point1 = '0 0 0'
#     point2 = '0 0 1'
#   []
#   [stress_pp]
#     type = RankTwoScalarAux
#     rank_two_tensor = stress
#     variable = stress_pp
#     scalar_type = HoopStress
#     point1 = '0 0 0'
#     point2 = '0 0 1'
#   []
# []

[Modules]
  [FluidProperties]
    [the_simple_fluid]
      type = SimpleFluidProperties
      bulk_modulus = 11E9
      viscosity = 1.0e-3
      density0 = 534
      thermal_expansion = 0.0000118
      cp = 4194
      cv = 4186
      porepressure_coefficient = 0
    []
  []
[]

# [ICs]
#   [tracer_concentration]
#     type = FunctionIC
#     # function = '0.5*(if x < 1.0001, 1, 0 )'
#     function = '0.001*if(x<1.0001,1,0)'
#     variable = tracer_concentration
#   []
# []

[Materials]
  [porosity]
    type = PorousFlowPorosity
    porosity_zero = 0.99
  []
  [biot_modulus]
    type = PorousFlowConstantBiotModulus
    solid_bulk_compliance = 2E-7
    fluid_bulk_modulus = 1E7
  []
  [permeability_aquifer]
    type = PorousFlowPermeabilityConst
    # block = aquifer
    permeability = '1E-14 0 0   0 1E-14 0   0 0 1E-14'
  []
  # [permeability_caps]
  #   type = PorousFlowPermeabilityConst
  #   block = caps
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
    block = 'caps aquifer'
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
  compute_scaling_once = true
  solve_type = Newton
  end_time = 50000
  dt = 100
  nl_abs_tol = 1E-10
  nl_rel_tol = 1E-14
  resid_vs_jac_scaling_param = 0.5
[]

[Outputs]
  exodus = true
  [out]
    type = Checkpoint
    interval = 5
    additional_execute_on = 'FINAL'
  []
[]
