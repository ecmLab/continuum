[Mesh]
  patch_size = 80
  patch_update_strategy = auto
  parallel_type = REPLICATED
  [./mesh]
    type = GeneratedMeshGenerator
    xmax = 5.0
    xmin = 0.0
    ymax = 1.0
    ymin = 0.0
    dim = 2
    elem_type = QUAD4
    nx = 5
    ny = 2
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
#    block = 'interLayer'
  [../]
[]

[AuxVariables]
  [flux_x]
    order = CONSTANT
    family = MONOMIAL
#    block = interLayer
  []
  [flux_y]
    order = CONSTANT
    family = MONOMIAL
#    block = interLayer
  []

[]

[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress hydrostatic_stress'
    use_automatic_differentiation = true
  [../]
[]

[AuxKernels]
  [./flux_x]
    type = AnisoTropicDiffusionFluxAux
    variable = flux_x
    diffusion_variable = conc
    component = x
  [../]
  [./flux_y]
    type = AnisoTropicDiffusionFluxAux
    variable = flux_y
    diffusion_variable = conc
    component = y
  [../]
[]

[Kernels]
  [./diffusion]
    type = ADMatAnisoDiffusion
    variable = conc
    diffusivity = diffusivity
    use_displaced_mesh = false
#    block = 'interLayer'
  [../]
  [./li_metal_dt]
    type = ADTimeDerivative
    variable = conc
    use_displaced_mesh = false
#    block = 'interLayer'
  [../]

[]

[Materials]
  [./diffusivity_Li1]
    type = ADDiffusionAlongPrincipalDirectionsMaterial
#   block = 'interLayer'
    diffusivity_vector = '1e8 0 0'
  [../]


  [./elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 5e-3
    poissons_ratio = 0.3
  [../]
  [./stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    # block = '1'
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
    # block = '1'
    omega = 1e13
    alpha = '1 0 0'
    concentration = conc
    cref = 0.0
    intBnd = 'bottom'
    # internal_solve_full_iteration_history = true
    # internal_solve_output_on = always
  [../]
[]

[BCs]
  [./top_y]
    type = ADDirichletBC
    boundary = 'top'
    variable = disp_y
    value = 0.0
    preset = true
  [../]

  [./top_x]
    type = ADDirichletBC
    boundary = 'top'
    variable = disp_x
    value = 0.0
    preset = true
  [../]
  [./top_flux]
    type = ADNeumannBC
    variable = conc
    value = 1e-14
    boundary = bottom
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
  petsc_options_iname = '-pc_type -pc_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'
  l_max_its = 50
  nl_max_its = 25

  nl_rel_tol = 1e-10

  start_time = 0.0
  dt = 0.001
  dtmax = 2.0
  dtmin = 1e-5
  num_steps = 10

[] # Executioner

[Outputs]
  [./out]
    type = Exodus
    file_base = rst/flat
  [../]
[] # Outputs

[Debug]
  show_var_residual_norms = true
  show_material_props = true
[]
