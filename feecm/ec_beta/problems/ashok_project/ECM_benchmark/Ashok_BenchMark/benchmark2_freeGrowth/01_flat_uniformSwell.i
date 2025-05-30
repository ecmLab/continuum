m1 = 0.5
m2 = 0.5
m3 = 0
out = arial
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
## Need to define this block otherwise simulation won't run
## And Ask DComputeFiniteStrain
[Problem]
  coord_type = RZ
[]
[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [./disp_x]
  [../]

  [./disp_y]
  [../]

[]
## AuxVariables are feeded into ADIsoTropicHyperViscoSwelling not
## boundry condtions
[AuxVariables]
  [./conc]
  [../]
[]
[AuxKernels]
  [./conc]
    type = FunctionAux
    variable = conc
    function = '1.0e-13*t'
    #block = 'interLayer'
  [../]
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
[Materials]
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
    alpha = '${fparse m1} ${fparse m2} ${fparse m3}'
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
    boundary = 'bottom'
    variable = disp_y
    value = 0.0
    preset = true
  [../]

  [./top_x]
    type = ADDirichletBC
    boundary = 'bottom'
    variable = disp_x
    value = 0.0
    preset = true
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
  dtmin = 1e-9
  num_steps = 1000

[] # Executioner

[Outputs]
  [./out]
    type = Exodus
    file_base = rst/01_flat_unifrom_visco_swell_${out}
  [../]
[] # Outputs

[Debug]
  show_var_residual_norms = true
  show_material_props = true
[]
