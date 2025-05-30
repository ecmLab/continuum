# Basic model to test directional swelling
# This file can be used for template for all types of swelling, can be simply
# done using the test system (TBD)
# Simple block model consiting of 2 blocks
# 1-> Interphase layer capable of directional swelling and plasticity
# 2-> Metal layer on top with plasticity
# Boundary conditions
#  left -> fixed in x-direction
#  right -> fixed in y-direction
# Swelling is induced by a linear concentration variation in time
# Checks to be performed
#   1) Fiber direction swelling with assumed principal directions
#   2) Areal swelling with assumed principal directions
#   3) Isotropic swelling with assumed principal directions
# All of these tests should be repeated by supplying the boundary names (intBnd)
# in the swelling block to auto-calculate the normals and principal directions
m1 = 0.5
m2 = 0.5
m3 = 0
out  = areal
[Mesh]
  [./mesh]
    type = GeneratedMeshGenerator
    nx = 10
    ny  = 5
    xmax = 5
    xmin = 0
    ymin = 0
    ymax = 0.5
    dim = 2
  [../]
[]

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


[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    #use_displaced_mesh = true
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy vonmises_stress strain_xx strain_yy elastic_strain_xx elastic_strain_yy'
    use_automatic_differentiation = true
  [../]
[]



[AuxVariables]
  [./conc]

  [../]

  [./Fp_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plastic_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strength]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./conc]
    type = FunctionAux
    variable = conc
    function = '8.32e-3*t'
  [../]

  [./fp_yy]
    type = ADRankTwoAux
    variable = Fp_yy
    rank_two_tensor = plastic_distortion
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./peeq]
    type = ADMaterialRealAux
    variable = plastic_strain
    property = effective_plastic_strain
    execute_on = timestep_end
  [../]
  [./strength]
    type = ADMaterialRealAux
    variable = strength
    property = yield_strength
    execute_on = 'TIMESTEP_END'
  [../]
[]

[BCs]
  [./symmy]
    type = ADDirichletBC
    variable = disp_x
    boundary = 'left right'
    value = 0
    preset = true
  [../]
  [./boty]
    type = ADDirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
    preset = true
  [../]
[]


[Materials]

  [./elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 7.83e3
    poissons_ratio = 0.3
  [../]
  [./stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
  [../]
  [./plas]
    type = ADIsotropicElasticSwelling
    absolute_tolerance = 1e-8
    relative_tolerance = 1e-08
    isotropic_swelling = true

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 1.0
    omega = 13
    concentration = conc
    cref = 0
    alpha = '${m1} ${m2} ${m3}'
    intBnd = 'bottom'
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
  [../]
[]

[Postprocessors]
  [./stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  [../]
  [./strain_yy]
    type = ElementAverageValue
    variable = strain_yy
  [../]
  # [./u_y]
  #   type = AverageNodalVariableValue
  #   variable = disp_y
  #   boundary = top
  # [../]
  [./peeq]
    type = ElementAverageValue
    variable = plastic_strain
  [../]
  [./fp_yy]
    type = ElementAverageValue
    variable = Fp_yy
  [../]
  [./stregnth]
    type = ElementAverageValue
    variable = strength
  [../]
  [./elastic_strain_yy]
    type = ElementAverageValue
    variable = elastic_strain_yy
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient

  #Preconditioned JFNK (default)
  automatic_scaling = true
  dt = 0.5
  solve_type = 'NEWTON'
#  petsc_options_iname = '-pc_type -pc_mat_solver_package'
#  petsc_options_value = 'lu superlu_dist'
  petsc_options = '-snes_converged_reason -ksp_converged_reason'
  petsc_options_iname = '-pc_type -mat_mffd_err -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu       1e-5          NONZERO               1e-15'
  dtmax = 5
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-11
  nl_max_its = 100
  dtmin = 1.0e-5
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.5
    growth_factor = 1.2
    cutback_factor = 0.65
    optimal_iterations = 25
  [../]
  # num_steps = 10
  end_time = 50.0
[]

[Outputs]
  # exodus = true
  # csv = true
  [./out]
    type = Exodus
    file_base = rst/Benchmark2_UG_elastic_swelling_${out}
  [../]
[]

[Debug]
  show_material_props = true
[]
