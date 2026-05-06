# Darcy flow with heat advection and conduction, and elasticity
# Dissolution of Ag in Li is the problem we are trying to solve here
# Model consists of a single block of porous material
#   1) Single phase fluid (Li/ LiAg)
#   2) Ag dissolves in the fluid
# If we want to model Ag precipitation in the fluid then perhaps we need 2 phases
# Li metal viscosity is set to 1e4
# Conversion from current density to mass flux = i A/m^2 / F(C/mol) * rho (kg/m^2) * omega (m^3/mol) = kg/m^2s
#                                               = current_density / Faradays' * molar volume * density
#                                               = 1 A/m^2 / 96485.2239 (C/mol) * 1.2e-5 (m^3/mol) * 534 (kg/m^3)
#                                               = 66.143144e-9 kg/m^2s
# Permeability in the Porous layer is Carman Kozeny formulation of type 1 where k_ij = A k_ij^0 * phi^n/ (1-phi)^m
#                                                                   A = f d^2 -> d is a grain size or particle size
# [Mesh]
#   [file]
#     type = FileMeshGenerator
#     file = test.e
#   []
#   [secondary_boundary_block]
#     type = LowerDBlockFromSidesetGenerator
#     input = file
#     sidesets = 'p_bottom'
#     new_block_name = 'secondary'
#   []
#   [primary_boundary_block]
#     type = LowerDBlockFromSidesetGenerator
#     input = secondary_boundary_block
#     sidesets = 'e_top'
#     new_block_name = 'primary'
#   []
# []
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
    boundary_name_prefix = 'p'
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
    new_block_id = '100 101'
    input = 'outflow_area'
  []
  [rename2]
    type = RenameBlockGenerator
    old_block_id = '100 101'
    input = rename
    new_block_name = 'LiLayer AgC'
  []

  [electrolyte]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 45
    ny = 15
    xmin = 0.0
    xmax = 100e-6
    ymin = -20e-6
    ymax = 0.0
    boundary_name_prefix = 'e'
    boundary_id_offset = 100
  []
  [rename_electrolyte]
    type = RenameBlockGenerator
    old_block_id = 0
    new_block_id = 200
    input = electrolyte
  []
  [full_mesh]
    type = MeshCollectionGenerator
    inputs = 'rename2 rename_electrolyte'
  []
  # [full_mesh]
  #   type = StitchedMeshGenerator
  #   inputs = 'rename2 rename_electrolyte'
  #   stitch_boundaries_pairs = 'p_bottom e_top'
  #   clear_stitched_boundary_ids = false
  # []
  [rename_electrolyte2]
    type = RenameBlockGenerator
    old_block_id = 200
    new_block_name = 'electrolyte'
    input = full_mesh
  []
  [secondary_boundary_block]
    type = LowerDBlockFromSidesetGenerator
    input = rename_electrolyte2
    sidesets = 'p_bottom'
    new_block_name = 'secondary'
  []
  [primary_boundary_block]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_boundary_block
    sidesets = 'e_top'
    new_block_name = 'primary'
  []
  # [interface]
  #   type = BreakMeshByBlockGenerator
  #   input = rename_electrolyte2
  #   block_pairs = '101 200'
  #   add_interface_on_two_sides = true
  #   split_interface = true
  #   show_info = true
  # []
  # [common_boundary]
  #   type = RenameBoundaryGenerator
  #   input = rename_electrolyte2
  #   old_boundary = 'p_bottom e_top'
  #   new_boundary = 'interface interface'
  #   show_info = true
  # []

[]

[Problem]
  coord_type = RZ
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  acceptable_iterations = 2
  group_variables = 'disp_r disp_z normal_lm; li_ion_V thermal_lm'

[]

[GlobalParams]
  PorousFlowDictator = dictator
  biot_coefficient = 1.0
  displacements = 'disp_r disp_z'
  gravity = '0 0 0'
[]

[Variables]
  [disp_r]
  []
  [disp_z]
  []
  [normal_lm]
    block = 'secondary'
  []
  [thermal_lm]
    block = 'secondary'
  []
  [li_ion_V]
    block = 'AgC LiLayer electrolyte'
  []
[]

[Constraints]
  [thermal_constraint]
    type = GapDisplacementConductanceConstraint
    variable = thermal_lm
    secondary_variable = li_ion_V
    primary_boundary = 'e_top'
    secondary_subdomain = 'secondary'
    secondary_boundary = 'p_bottom'
    primary_subdomain = 'primary'
    # k_function = gapk
    k = 1e2
    use_displaced_mesh = true
    compute_lm_residuals = true
    extra_vector_tags = 'ref'
    include_gap = false
    include_equilibrium_potential = false
    R = 8.3145
    faraday = 96485.3329
    temperature = 298
    surfaceType = SECONDARY
    include_concentration = false

  []
  [normal_lm]
    type = NormalNodalLMMechanicalContact
    secondary = 'p_bottom'
    primary = 'e_top'
    variable = normal_lm
    primary_variable = disp_r
    normal_smoothing_distance = 0.2
    disp_y = disp_z
    ncp_function_type = min
    tangential_tolerance = 0.2
    # c = 1e-3
  []
  [normal_r]
    type = NormalMortarMechanicalContact
    primary_boundary = 'e_top'
    secondary_boundary = 'p_bottom'
    primary_subdomain = 'primary'
    secondary_subdomain = 'secondary'
    variable = normal_lm
    secondary_variable = disp_r
    component = x
    use_displaced_mesh = true
    compute_lm_residuals = false
    # extra_vector_tags = 'ref'
  []
  [normal_z]
    type = NormalMortarMechanicalContact
    secondary_boundary = 'e_top'
    primary_boundary = 'p_bottom'
    primary_subdomain = 'primary'
    secondary_subdomain = 'secondary'
    variable = normal_lm
    secondary_variable = disp_z
    component = y
    use_displaced_mesh = true
    compute_lm_residuals = false
    # extra_vector_tags = 'ref'
  []
[]

[Modules/TensorMechanics/Master]
  [all]
    add_variables = true
    strain = SMALL
    volumetric_locking_correction = false
    generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress hydrostatic_stress'
    use_automatic_differentiation = false
    # block = 'AgC LiLayer electrolyte'
    # block = 'AgC LiLayer electrolyte'
    extra_vector_tags = 'ref'
    # use_finite_deform_jacobian = true
  []
[]

[Kernels]
  [li_ion_V]
    type = ADHeatConduction
    block = 'AgC LiLayer electrolyte'
    variable = li_ion_V
    use_displaced_mesh = false
  []
[]

[BCs]
  [constant_current_density]
    type = ADFunctionNeumannBC
    boundary = 'e_bottom'
    variable = li_ion_V
    function = 'if (t <= 10, 0.1*t, 1.0 )'
    extra_vector_tags = 'ref'
  []
  [constant_Voltage]
    type = DirichletBC
    boundary = 'p_top'
    variable = li_ion_V
    value = 0.0
  []
  [pinned_top_bottom_r]
    type = DirichletBC
    variable = disp_r
    value = 0
    boundary = 'p_left e_left' # p_right e_right'
  []
  [pinned_top_bottom_z]
    type = DirichletBC
    variable = disp_z
    value = 0
    boundary = 'p_top'
  []
  [pressure_bottom]
    type = ADPressure
    variable = disp_z
    component = 1
    extra_vector_tags = 'ref'
    function = '1e6'
    # function = 'if (t <=10, 0.1e5*t, 1e5)'
    boundary = 'e_bottom'
  []
[]

[Materials]
  [thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1.0e-2
    # block = 'blockCeramic'
  []

  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 5E9
    poissons_ratio = 0.25
    # block = 'AgC LiLayer electrolyte'
  []
  [stress]
    type = ComputeLinearElasticStress
    # block = 'AgC LiLayer electrolyte'
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
    # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_force_iteration'
    # petsc_options_value = ' lu       mumps 1'
  []
[]

[Executioner]
  type = Transient
  automatic_scaling = true
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor'
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration '
                        '-pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1     NONZERO               1e-20               '

  # compute_scaling_once = false
  solve_type = Newton
  # end_time = 36000
  dt = 10
  nl_max_its = 100
  nl_abs_tol = 1E-8
  nl_rel_tol = 1E-5
  resid_vs_jac_scaling_param = 0.5
  dtmax = 10
  end_time = 500
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
  execute_on = 'INITIAL TIMESTEP_END'
  # sync_times = '500'
  # sync_only = false
  file_base = orig
[]
