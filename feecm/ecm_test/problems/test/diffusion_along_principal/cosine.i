[Mesh]
  patch_size = 80
  patch_update_strategy = auto
  parallel_type = REPLICATED
  [./mesh]
    type = FileMeshGenerator
    file = cosine_conformal.msh
  [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [./disp_x]
  [../]

  [./disp_y]
  [../]

  [./conc]
    initial_condition = 0
    block = 'interLayer'
  [../]
[]
# [AuxKernels]
#   [./flux_x]
#     type = AnisoTropicDiffusionFluxAux
#     diffusivity = diffusivity
#     variable = flux_x
#     component = 0
#   [../]
#   [./flux_x]
#     type = AnisoTropicDiffusionFluxAux
#     diffusivity = diffusivity
#     variable = flux_y
#     component = 1
#   [../]
# []

[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress hydrostatic_stress'
    use_automatic_differentiation = true
  [../]
[]

[Kernels]
  [./diffusion]
    type = ADMatAnisoDiffusion
    variable = conc
    diffusivity = diffusivity
    use_displaced_mesh = false
    block = 'interLayer'
  [../]
  [./li_metal_dt]
    type = ADTimeDerivative
    variable = conc
    use_displaced_mesh = false
    block = 'interLayer'
  [../]

[]

[Materials]
  [./diffusivity_Li1]
    type = ADDiffusionAlongPrincipalDirections
    diffusivity_name =  'diffusivity'
    diffusivity_vector = '1e8 0 0'
    block = 'interLayer'
    # output_properties = diffusion_tensor
  [../]

  [./elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 5e-3
    poissons_ratio = 0.3
  [../]
  [./stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    block = 'interLayer'
  [../]
  [./plas]
    type = ADIsoTropicHyperViscoSwelling
    # absolute_tolerance = 1e-5
    hardening_exponent = 1.8
    saturation_resistance = 8.0e-6
    initial_resistance = 2.0e-6
    hardening_modulus = 40.0e-6
    rate_exponent = 0.18
    reference_strain_rate = 0.05
    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 1.0
    omega = 1e13
    alpha = '1 0 0'
    concentration = conc
    cref = 0.0
    intBnd = 'blockMetal_bottom'
    block = 'interLayer'
    # internal_solve_full_iteration_history = true
    # internal_solve_output_on = always
  [../]
  [./stress2]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas2'
    block = 'blockMetal'
  [../]
  [./plas2]
    type = ADIsoTropicHyperVisco
    # absolute_tolerance = 1e-5
    hardening_exponent = 1.8
    saturation_resistance = 8.0e-6
    initial_resistance = 2.0e-6
    hardening_modulus = 40.0e-6
    rate_exponent = 0.18
    reference_strain_rate = 0.05
    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 1.0
    block = 'blockMetal'
  [../]

[]

[BCs]
  [./top_y]
    type = ADDirichletBC
    boundary = 'blockMetal_top'
    variable = disp_y
    value = 0.0
    preset = true
  [../]

  [./top_x]
    type = ADDirichletBC
    boundary = 'blockMetal_top'
    variable = disp_x
    value = 0.0
    preset = true
  [../]
  [./top_flux]
    type = ADNeumannBC
    variable = conc
    value = 1e-14
    boundary = 'blockMetal_bottom'
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
  solve_type = 'NEWTON'
  automatic_scaling = true
  compute_scaling_once = false
  petsc_options = '-snes_converged_reason -ksp_converged_reason'
  # petsc_options_iname = '-pc_type -pc_mat_solver_package'
  # petsc_options_value = 'lu superlu_dist'
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1     NONZERO               1e-20               '

  l_max_its = 50
  nl_max_its = 25
  l_tol = 1e-03
  nl_abs_tol = 1e-15
  nl_rel_tol = 1e-6

  start_time = 0.0
  dt = 1.0e-4
  dtmax = 2.0
  dtmin = 1e-5
  # num_steps = 10
  end_time = 10.0
 scaling_group_variables = 'disp_x disp_y conc'
 [./TimeStepper]
   type = IterationAdaptiveDT
   dt = 1e-4
   growth_factor = 1.5
   cutback_factor = 0.5
   optimal_iterations = 40
   # timestep_limiting_postprocessor = matl_ts_min
 [../]

[] # Executioner

[Outputs]
  [./out]
    type = Exodus
    # output_material_properties = true
  [../]
[] # Outputs

[Debug]
  show_var_residual_norms = true
  show_material_props = true
[]
