
[Mesh]
  [./mesh1]
    type = GeneratedMeshGenerator
    xmax = 1e-3
    xmin = -1.0e-3
    ymax = 0.5e-3
    ymin = 0.0
    dim = 2
    nx = 25
    ny = 15
  [../]
  [./m1]
    type = RenameBlockGenerator
    input = mesh1
    old_block_id = '0'
    new_block_id = '1'
  [../]
  [./slave_boundary]
    type = RenameBoundaryGenerator
    input = m1
    old_boundary_name = 'bottom right top left'
    new_boundary_id = '16 17 18 19'
  [../]

  [./mesh2]
    type = GeneratedMeshGenerator
    xmax = 1e-3
    xmin = -1.0e-3
    ymax = 0.0
    ymin = -0.5e-3
    dim = 2
    nx = 15
    ny = 15
  [../]

  [./m2]
    type = RenameBlockGenerator
    input = mesh2
    old_block_id = '0'
    new_block_id = '2'
  [../]

  [./master_boundary]
    type = RenameBoundaryGenerator
    input = m2
    old_boundary_name = 'bottom right top left'
    new_boundary_id = '11 12 15 14'
  [../]
  [./cmg]
    type = CombinerGenerator
    inputs = 'master_boundary slave_boundary'
  [../]
  [./slave_boundary_block]
    type = LowerDBlockFromSidesetGenerator
    input = cmg
    sidesets = '16'
    new_block_id = '4'
  [../]
  [./master_boundary_block]
    type = LowerDBlockFromSidesetGenerator
    input = slave_boundary_block
    sidesets = '15'
    new_block_id = '3'
  [../]


[]


[GlobalParams]
  displacements = 'ux uy'
[]

[Variables]
  [./ux]
    block = '1 2'
  [../]
  [./uy]
    block = '1 2'
  [../]
  [./temp]
    block = '1 2'
  [../]
  [./thermal_lm]
    block = '4'
  [../]
[]

# [ICs]
#   [./uy]
#     type = ConstantIC
#     variable = uy
#     value = -1.0e-3
#     block = 1
#   [../]
# []

# [AuxVariables]
#   [./conc]
#     block = '2'
#   [../]
# []
# [AuxKernels]
#   [./conc]
#     type = FunctionAux
#     variable = conc
#     function = '1.0e-15*t'
#     block = '2'
#   [../]
# []
#
[Functions]
  [./pressure]
    type = ConstantFunction
    value = 1.0e-6
  [../]
  [./uy]
    type = ParsedFunction
    value = '-1.0e-3*t'
  [../]
[]

# [Modules/TensorMechanics/Master]
#   [./all]
#     add_variables = true
#     strain = FINITE
#     use_displaced_mesh = true
#     # volumetric_locking_correction = true
#     generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress'
#     use_automatic_differentiation = true
#     block = '1 2'
#     # use_finite_deform_jacobian = true
#   [../]
# []
# [Contact]
#   [./mech_contact]
#     disp_x = ux
#     disp_y = uy
#     master = 15
#     slave = 16
#     penalty = 1e-2
#   [../]
# []
[Constraints]
  [./thermal_constraint]
    type = GapDisplacementConductanceConstraint
    variable = thermal_lm
    slave_variable = temp
    slave_boundary = 16
    slave_subdomain = 4
    master_boundary = 15
    master_subdomain = 3
    gap_limit = 1.0e-9
    k = 1000.0
    # use_displaced_mesh = true
  [../]
[]
[Kernels]
  [./temp]
    type = ADHeatConduction
    diffusion_coefficient = thermal_conductivity
    block = '1 2'
    variable = temp
  [../]
  # [./temp_dt]
  #   type = ADHeatConductionTimeDerivative
  #   density_name = density
  #   variable = temp
  #   block = '1 2'
  #   specific_heat = 1.0
  # [../]
[]
[Materials]
  [./thermal_conductivity]
    type = HeatConductionMaterial
    thermal_conductivity = 5.0
    block = '1 2'
  [../]
  [./density]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = '1.0'
    block = '1 2'
  [../]

  # [./elasticity_tensor_Li]
  #   type = ComputeIsotropicElasticityTensor
  #   youngs_modulus = 7.81e-3
  #   poissons_ratio = 0.38
  #   block = '2'
  # [../]
  # [./elasticity_tensor_Llzo]
  #   type = ComputeIsotropicElasticityTensor
  #   youngs_modulus = 170e-3
  #   poissons_ratio = 0.3
  #   block = '1'
  # [../]
  # [./stress_llzo]
  #   type = ADComputeFiniteStrainElasticStress
  #   block = '1'
  # [../]
  #
  # [./stress_Li]
  #   type = ADComputeFiniteStrainElasticStress
  #   # inelastic_models = 'plas'
  #   # perform_finite_strain_rotations = true
  #   block = '2'
  # [../]
  # [./plas]
  #   type = ADC
  # [../]
[]

[BCs]
  # [./Li_top]
  #   type = ADFunctionPresetBC
  #   boundary = '18'
  #   variable = uy
  #   component = 1
  #   function = uy
  # [../]
  #
  # [./left_right]
  #   type = ADPresetBC
  #   variable = ux
  #   boundary = '11 12 14 17 19'
  #   value = 0
  # [../]
  # [./bottom]
  #   type = ADPresetBC
  #   variable = uy
  #   boundary = '11'
  #   value = 0
  # [../]
  [./temp_top]
    type = ADPresetBC
    boundary = '18'
    variable = temp
    value = 100.0
  [../]
  [./temp_bot]
    type = ADPresetBC
    boundary = '11'
    variable = temp
    value = 10.0
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./left]
    type = SideFluxIntegral
    variable = temp
    boundary = 15
    diffusivity = thermal_conductivity
  [../]
  [./right]
    type = SideFluxIntegral
    variable = temp
    boundary = 16
    diffusivity = thermal_conductivity
  [../]
[]

[Executioner]

  type = Transient

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'
  automatic_scaling = true
  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre    boomeramg      101'
  # petsc_options_iname = -pc_type
  # petsc_options_value = lu


  line_search = 'none'


  nl_abs_tol = 5e-8
  # nl_rel_tol = 1e-6

  l_max_its = 100


  start_time = 0.0
  dt = 0.001
  dtmax = 0.05
  dtmin = 1e-5
  end_time = 2.0
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.001
    growth_factor = 1.5
    cutback_factor = 0.25
    optimal_iterations = 40
    # timestep_limiting_postprocessor = matl_ts_min
  [../]

[] # Executioner

[Outputs]
  [./out]
    type = Exodus
    elemental_as_nodal = true
  [../]
  [./csv]
    type = CSV
  [../]
[] # Outputs
