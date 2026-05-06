# This test involves only thermal expansion strains on a 2x2x2 cube of approximate
# steel material; however, in this case the stress free temperature of the material
# has been set to 200K so that there is an initial delta temperature of 100K.

# An initial temperature of 0K is given for the material,
# and an auxkernel is used to calculate the temperature in the entire cube to
# raise the temperature each time step.  The final temperature is 90K
# The thermal strain increment should therefore be
#     (90K - 0K) * 7.5e-4 1/K = 6.75e-2 m/m.

# This test uses a start up step to identify problems in the calculation of
# eigenstrains with a stress free temperature that is different from the initial
# value of the temperature in the problem
#by: Shafee Farzanian
# units are: Length=um, Force=uN, Mass=?, Time=hour, Stress=MPa, Energy=?

[Mesh]
  patch_size = 80
  patch_update_strategy = auto
  [./fmg]
    type = FileMeshGenerator
    file = data/noCoating.msh
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  #volumetric_locking_correction = true
[]


[AuxVariables]
  [./temp]
    initial_condition = 0
  [../]
  [./eigenstrain_yy]
    order = FIRST
    family = MONOMIAL
    block = 'block_NMC '
  [../]
  [./eigenstrain_xx]
    order = FIRST
    family = MONOMIAL
    block = 'block_NMC '
  [../]
  [./total_strain_yy]
    order = FIRST
    family = MONOMIAL
    block = 'block_NMC '
  [../]
  [./total_strain_xx]
    order = FIRST
    family = MONOMIAL
    block = 'block_NMC '
  [../]
[]

[Functions]
  [./temperature_load]
    type = ParsedFunction
    #value = t*180
    value = 'if(t<=0.5, 180*t, (-1)*(t-0.5)*180 + 90.0)'
  [../]

  [./hf]
    type = PiecewiseLinear
    x = '0  '
    y = '12600 '
  [../]

  [./hf1]
    type = PiecewiseLinear
    x = '0  0.1'
    y = '500  500'
  [../]

[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./NMC]
        strain = FINITE
        #use_displaced_mesh = true
        #incremental = true
        add_variables = true
        eigenstrain_names = eigenstrain
        generate_output = 'stress_xx stress_yy vonmises_stress strain_xx strain_yy'
        block = 'block_NMC '
        #extra_vector_tags = 'ref'
        #use_automatic_differentiation = true
      [../]
      [./LPS]
      strain = FINITE
      #use_displaced_mesh = true
      add_variables = true
      generate_output = 'stress_yy stress_xx vonmises_stress plastic_strain_xx plastic_strain_yy'
      block = ' block_LPS'

      [../]
    [../]
  [../]
[]

[AuxKernels]
  [./tempfuncaux]
    type = FunctionAux
    variable = temp
    function = temperature_load
  [../]

  [./eigenstrain_yy]
    type = RankTwoAux
    block = 'block_NMC '
    rank_two_tensor = eigenstrain
    variable = eigenstrain_yy
    index_i = 1
    index_j = 1
    execute_on = 'initial timestep_end'
  [../]
  [./eigenstrain_xx]
    type = RankTwoAux
    block = 'block_NMC '
    rank_two_tensor = eigenstrain
    variable = eigenstrain_xx
    index_i = 0
    index_j = 0
    execute_on = 'initial timestep_end'
  [../]

  [./total_strain_yy]
    type = RankTwoAux
    block = 'block_NMC '
    rank_two_tensor = total_strain
    variable = total_strain_yy
    index_i = 1
    index_j = 1
    execute_on = 'initial timestep_end'
  [../]
  [./total_strain_xx]
    type = RankTwoAux
    block = 'block_NMC '
    rank_two_tensor = total_strain
    variable = total_strain_xx
    index_i = 0
    index_j = 0
    execute_on = 'initial timestep_end'
  [../]
[]

[Contact]
  [nmc_lps]
    primary = 'block_NMC_right'
    secondary = 'block_LPS_left'

    penalty = 1e5
    formulation = penalty
    tangential_tolerance = 0.0001

  []
[]

[BCs]
  [./x_disp]
    type = DirichletBC
    variable = disp_x
    boundary = 'block_right block_top block_left'
    value = 0.0
  [../]
  [./y_disp]
    type = DirichletBC
    variable = disp_y
    boundary = 'block_right block_top block_bottom'
    value = 0.0
  [../]

[]

[Materials]

  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    block = 'block_NMC'
    youngs_modulus = 177.5e3
    poissons_ratio = 0.33
  [../]
  [./small_stress]
    type = ComputeFiniteStrainElasticStress
    block = 'block_NMC'
  [../]
  [./thermal_expansion_strain]
    type = ComputeThermalExpansionEigenstrain
    block = 'block_NMC'
    stress_free_temperature = 0
    thermal_expansion_coeff = 0.00075
    temperature = temp
    eigenstrain_name = eigenstrain
  [../]

  [./elasticity_tensor1]
    type = ComputeIsotropicElasticityTensor
    block = 'block_LPS'
    youngs_modulus = 18.5e3
    poissons_ratio = 0.33
    #eigenstrain_name = eigenstrain
  [../]
  #[./small_stress1]
  [./isotropic_plasticity1]
    type = IsotropicPlasticityStressUpdate
    block = 'block_LPS'
    yield_stress = 500.0
    hardening_function = hf1
  [../]
  [./radial_return_stress1]
    type = ComputeMultipleInelasticStress
    tangent_operator = elastic
    inelastic_models = 'isotropic_plasticity1'
    block = 'block_LPS'
  [../]

[]

[Executioner]
  type = Transient
  automatic_scaling = true
  solve_type = NEWTON
  petsc_options_iname = -pc_type
  petsc_options_value = lu
#  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart'
#  petsc_options_value = 'asm lu 1 101'

  nl_max_its = 99
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-9
  l_tol = 1e-8

  start_time = 0.0
  n_startup_steps = 1
  end_time = 1
  dt = 0.025
  dtmin = 0.001
[]

[Outputs]
  exodus = true
  file_base = rst/rst
[]

[Postprocessors]
  [./eigenstrain_xx]
    type = ElementAverageValue
    variable = eigenstrain_xx
    execute_on = 'initial timestep_end'
    block = 'block_NMC '

  [../]
  [./eigenstrain_yy]
    type = ElementAverageValue
    variable = eigenstrain_yy
    execute_on = 'initial timestep_end'
    block = 'block_NMC '

  [../]
  [./total_strain_xx]
    type = ElementAverageValue
    variable = total_strain_xx
    execute_on = 'initial timestep_end'
    block = 'block_NMC '

  [../]
  [./total_strain_yy]
    type = ElementAverageValue
    variable = total_strain_yy
    execute_on = 'initial timestep_end'
    block = 'block_NMC '

  [../]
  [./temperature]
    type = AverageNodalVariableValue
    variable = temp
    execute_on = 'initial timestep_end'
    block = 'block_NMC '

  [../]
[]
