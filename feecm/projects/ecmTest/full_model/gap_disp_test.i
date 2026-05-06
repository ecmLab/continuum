elem = QUAD4
order = FIRST
[Mesh]
  patch_size = 80
  patch_update_strategy = auto
  parallel_type = REPLICATED
  [./mesh]
    type = FileMeshGenerator
    file = test_separated.msh
  [../]
  [./secondary_boundary_block]
    type = LowerDBlockFromSidesetGenerator
    input = mesh
    sidesets = 'blockMetal_bottom'
    new_block_name = 'secondary'
  [../]
  [./primary_boundary_block]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_boundary_block
    sidesets = 'blockCeramic_top'
    new_block_name = 'primary'
  [../]
[]

[AuxVariables]
  [./ux]
    block = 'blockCeramic blockMetal'
  [../]
  [./uy]
    block = 'blockCeramic blockMetal'
    # iniital_condition = 1e-3
  [../]
  [./bndliflux]
    order = CONSTANT
    family = MONOMIAL
    block = 'blockCeramic blockMetal'
  [../]

[]

[AuxKernels]
  [./ux]
    type = ConstantAux
    value = 0
    variable = ux
  [../]
  [./uy]
    type = ConstantAux
    value = 0.0
    variable = uy
    block = 'blockMetal'
  [../]
  [./bnd_li_flux]
    type = DiffusionFluxNormalToBoundaryAux
    variable = bndliflux
    boundary = 'blockMetal_bottom blockCeramic_top'
    diffusion_variable = temp
    diffusivity = thermal_conductivity
  [../]

[]

[Variables]
  [./temp]
    block = 'blockCeramic blockMetal'
    initial_condition = 0
  [../]
  [./thermal_lm]
    block = 'secondary'
  [../]
[]
[Functions]
  [./gapk]
    type = PiecewiseLinear
    data_file = 'gap_cond.csv'
    # x = '-1.0 0.0 1e-6 2e-6 3e-6 4e-6 5e-6 1e-3'
    # y = '100.0 100.0 100.0 100.0 0.0 0.0 0.0 0.0'
    format = columns
  [../]
[]


[Constraints]
  [./thermal_constraint]
    type = GapDisplacementConductanceConstraint
    variable = thermal_lm
    secondary_variable = temp
    primary_boundary =  'blockCeramic_top'
    secondary_subdomain = 'secondary'
    secondary_boundary = 'blockMetal_bottom'
    primary_subdomain = 'primary'
    k_function = gapk
    k = 1e2
    use_displaced_mesh = true
    compute_lm_residuals = true
    displacements = 'ux uy'
  [../]
[]

[Kernels]
  [./temp_stead]
    type = ADHeatConduction
    block = 'blockCeramic blockMetal'
    variable = temp
    use_displaced_mesh = false
  [../]
  [./temp_dt]
    type = ADTimeDerivative
    block = 'blockCeramic blockMetal'
    variable = temp
    use_displaced_mesh = false
  [../]
[]

[Materials]
  [./thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1.0e2
    block = 'blockCeramic'
  [../]

  [./thermal_conductivity2]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1e2
    block = 'blockMetal'
    use_displaced_mesh = true
  [../]
[]

[BCs]
  # [./top_temp]
  #   type = ADDirichletBC
  #   variable = temp
  #   boundary = 'blockMetal_top'
  #   value = 0.0
  # [../]

  [./bot_temp]
    type = ADNeumannBC
    boundary = 'blockCeramic_bottom'
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

[Executioner]

  type = Transient

  #Preconditioned JFNK (default)
  solve_type = 'NEWTON'
  automatic_scaling = true
  compute_scaling_once = false
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor'
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1     NONZERO               1e-20               '
  line_search = contact
  l_max_its = 50
  nl_max_its = 25
  nl_abs_tol = 1e-13
  # nl_rel_tol = 1e-1
  start_time = 0.0
  dt = 1.0
  dtmax = 1.0
  dtmin = 1e-5
  end_time = 10.0

[] # Executioner

[Outputs]
  [./out]
    type = Exodus
    [../]
[] # Outputs

[Debug]
  show_var_residual_norms = true
  show_material_props = true
[]
