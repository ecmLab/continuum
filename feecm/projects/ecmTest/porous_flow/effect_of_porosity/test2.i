[Mesh]
  [mesh]
    # type = FileMeshGenerator
    # file = test_mesh.e
    type = GeneratedMeshGenerator
    dim = 2
    nx = 60
    xmin = 0.0
    xmax = 100e-6
    bias_x = 1
    ny = 30
    ymin = 0
    ymax = 30e-6
  []
  [AgC]
    type = SubdomainBoundingBoxGenerator
    block_id = 1
    bottom_left = '0.0 0.0 0'
    top_right = '100e-6 15.01e-6 0'
    input = mesh
  []
  [injection_area]
    type = ParsedGenerateSideset
    # combinatorial_geometry = 'x*x+y*y<1.01'
    combinatorial_geometry = 'y<0.0001e-6'
    included_subdomain_ids = 1
    new_sideset_name = 'injection_area'
    input = 'AgC'
  []
  [outflow_area]
    type = ParsedGenerateSideset
    combinatorial_geometry = 'y>29.99e-6'
    included_subdomain_ids = 1
    new_sideset_name = 'outflow_area'
    input = injection_area
  []
  [rename]
    type = RenameBlockGenerator
    old_block_id = '0 1'
    new_block_name = 'LiLayer AgC'
    input = 'outflow_area'
  []
  [interface]
    type = SideSetsBetweenSubdomainsGenerator
    new_boundary = 'interface'
    primary_block = 'AgC'
    paired_block = 'LiLayer'
    input = rename
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
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  acceptable_iterations = 2
  group_variables = 'disp_r disp_z temperature porepressure'
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
  residual_saturation = 0.001
  van_genuchten_alpha = 1e-6
  van_genuchten_m = 0.6
[]

[BCs]
  [constant_injection_porepressure]
    type = DirichletBC
    variable = porepressure
    value = 1e5
    boundary = top
  []

  [constant_injection_flux]
    type = PorousFlowSink
    variable = porepressure
    flux_function = "if (t <= 500,if (t <10, -66.41344e-9*3*t,-66.41344e-9*30), 66.41344e-9*30)"
    # flux_function = -664.13e-9
    boundary = bottom
    fluid_phase = 0
    use_relperm = true
    extra_vector_tags = 'ref'
  []
  [constant_injection_temperature]
    type = DirichletBC
    variable = temperature
    value = 293
    boundary = bottom
  []
  [left]
    type = DirichletBC
    variable = disp_r
    value = 0
    boundary = 'left'
  []
  [roller_top_bottom]
    type = DirichletBC
    variable = disp_z
    value = 0
    boundary = 'top bottom'
  []
  [cavity_pressure_x]
    type = Pressure
    boundary = bottom
    variable = disp_r
    component = 0
    factor = 1E6
    use_displaced_mesh = false
    extra_vector_tags = 'ref'
  []
  [cavity_pressure_y]
    type = Pressure
    boundary = bottom
    variable = disp_z
    component = 1
    factor = 1E6
    use_displaced_mesh = false
    extra_vector_tags = 'ref'
  []
[]

[Postprocessors]
  [li_mass]
    type = PorousFlowFluidMass
    PorousFlowDictator = dictator
    fluid_component = 0
  []
  [porepressure_bottom]
    type = SideAverageValue
    variable = porepressure
    boundary = bottom
  []
  [stress_zz_bottom]
    type = SideAverageValue
    variable = stress_yy
    boundary = bottom
  []
  [saturation]
    type = SideAverageValue
    variable = saturation0
    boundary = bottom
  []
  [total_current]
    type = PorousFlowFluidMass
    PorousFlowDictator = dictator
    fluid_component = 0
  []
  # [avg_permeability]
  #   type = ElementAverageValue
  #   variable = permeability
  #   block = AgC
  # []
[]

[Modules]
  [FluidProperties]
    [the_simple_fluid]
      type = SimpleFluidProperties
      bulk_modulus = 11E9
      viscosity = 1.0e2
      density0 = 534
      thermal_expansion = 0.0000118
      cp = 4194
      cv = 4186
      porepressure_coefficient = 0
    []
  []
[]

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
    porosity = 0.95
    chemical = false
    block = LiLayer
    # initial_mineral_concentrations = initial_and_reference_conc
    # reference_chemistry = initial_and_reference_conc
  []

  [biot_modulus]
    type = PorousFlowConstantBiotModulus
    solid_bulk_compliance = 2E-9
    fluid_bulk_modulus = 11E9
  []
  [permeability_AgC]
    #   type = PorousFlowPermeabilityConst
    #   permeability = '1E-25 0 0   0 1E-15 0   0 0 1E-25'
    #   block = 'AgC'
    # []
    type = PorousFlowPermeabilityKozenyCarman
    block = AgC
    k_anisotropy = '1E-25 0 0   0 1 0   0 0 1E-25'
    f = 0.01
    d = 50e-9
    poroperm_function = kozeny_carman_fd2
    m = 3
    n = 2
  []
  [permeability_LiLayer]
    type = PorousFlowPermeabilityKozenyCarman
    block = LiLayer
    k_anisotropy = '1 0 0   0 1 0   0 0 1'
    f = 0.01
    d = 50e-9
    poroperm_function = kozeny_carman_fd2
    m = 3
    n = 2
  []

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
    thermal_expansion_coeff = 1e-19 # this is the linear thermal expansion coefficient
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
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_force_iteration'
    petsc_options_value = ' lu       mumps 1'
  []
[]

[Executioner]
  type = Transient
  automatic_scaling = true
  # compute_scaling_once = false
  solve_type = Newton
  # end_time = 36000
  dt = 10
  nl_max_its = 100
  nl_abs_tol = 1E-8
  nl_rel_tol = 1E-4
  resid_vs_jac_scaling_param = 0.5
  dtmax = 10
  end_time = 1000
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0
    optimal_iterations = 100
    cutback_factor = 0.5
    growth_factor = 1.1
  []

[]

[Outputs]
  exodus = true
  csv = true
  sync_times = '500'
  sync_only = false
  file_base = xxx2
[]
