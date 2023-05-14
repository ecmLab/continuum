# This test involves only thermal expansion strains on a 2x2x2 cube of approximate
# steel material; however, in this case the stress free temperature of the material
# has been set to 200K so that there is an initial delta temperature of 100K.

# An initial temperature of 300K is given for the material,
# and an auxkernel is used to calculate the temperature in the entire cube to
# raise the temperature each time step.  The final temperature is 675K
# The thermal strain increment should therefore be
#     (675K - 300K) * 1.3e-5 1/K + 100K * 1.3e-5 1/K = 6.175e-3 m/m.

# This test uses a start up step to identify problems in the calculation of
# eigenstrains with a stress free temperature that is different from the initial
# value of the temperature in the problem
#by: Shafee Farzanian
# units are: Length=mm, Force=N, Mass=tonne, Time=s, Stress=MPa, Energy=mJ

[Mesh]
  patch_size = 80
  patch_update_strategy = auto
  [./fmg]
    type = FileMeshGenerator
    file = RVE_NoCoating.msh
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
    block = 'block_NMC block_LPS'
  [../]
  [./eigenstrain_xx]
    order = FIRST
    family = MONOMIAL
    block = 'block_NMC block_LPS'
  [../]
  [./total_strain_yy]
    order = FIRST
    family = MONOMIAL
    block = 'block_NMC block_LPS'
  [../]
  [./total_strain_xx]
    order = FIRST
    family = MONOMIAL
    block = 'block_NMC block_LPS'
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
    x = '0  0.03'
    y = '12600 12610'
  [../]

  [./hf1]
    type = PiecewiseLinear
    x = '0  0.03'
    y = '2000  2005'
  [../]

[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./all]
        strain = FINITE
        #use_displaced_mesh = true
        incremental = true
        add_variables = true
        eigenstrain_names = eigenstrain
        generate_output = 'stress_xx stress_yy vonmises_stress strain_xx strain_yy plastic_strain_xx plastic_strain_yy'
        block = 'block_NMC block_LPS'
        #extra_vector_tags = 'ref'
        #use_automatic_differentiation = true
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
    block = 'block_NMC block_LPS'
    rank_two_tensor = eigenstrain
    variable = eigenstrain_yy
    index_i = 1
    index_j = 1
    execute_on = 'initial timestep_end'
  [../]
  [./eigenstrain_xx]
    type = RankTwoAux
    block = 'block_NMC block_LPS'
    rank_two_tensor = eigenstrain
    variable = eigenstrain_xx
    index_i = 0
    index_j = 0
    execute_on = 'initial timestep_end'
  [../]

  [./total_strain_yy]
    type = RankTwoAux
    block = 'block_NMC block_LPS'
    rank_two_tensor = total_strain
    variable = total_strain_yy
    index_i = 1
    index_j = 1
    execute_on = 'initial timestep_end'
  [../]
  [./total_strain_xx]
    type = RankTwoAux
    block = 'block_NMC block_LPS'
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

    #secondary = 2
    #primary = 3
    #model = frictionless
    #penalty = 1e+6
    #normalize_penalty = true
    #formulation = kinematic
    #normal_smoothing_distance = 0.1

    #primary = 20
    #secondary = 10
    #penalty = 1e8
    #formulation = penalty
    #tangential_tolerance = 1e-3
    #tension_release = -1

    #penalty = 1e5
    #formulation = kinematic

    #model = frictionless
    #formulation = penalty
    #penalty = 1e9
    #normalize_penalty = true

    #primary = 'bottom_top'
    #secondary = 'top_bottom'
    #formulation = mortar
    #model = coulomb
    #friction_coefficient = 0.4
    #c_normal = 1e4
    #c_tangential = 1.0e4
    #interpolate_normals = false

    #model = frictionless
    #formulation = mortar
    #c_normal = 1e-1

    penalty = 1e7
    formulation = penalty
    tangential_tolerance = 0.0001

    #model = frictionless
    #penalty = 1e+6
    #normal_smoothing_distance = 0.1
  []
[]

[BCs]
  [./x_left]
    type = DirichletBC
    variable = disp_x
    boundary = 'block_left'
    value = 0.0
  [../]
  [./y_bot]
    type = DirichletBC
    variable = disp_y
    boundary = 'block_bottom'
    value = 0.0
  [../]
  [./x_right]
    type = DirichletBC
    variable = disp_x
    boundary = 'block_right'
    value = 0.0
  [../]
  [./y_right]
    type = DirichletBC
    variable = disp_y
    boundary = 'block_right'
    value = 0.0
  [../]
  [./x_top]
    type = DirichletBC
    variable = disp_x
    boundary = 'block_top'
    value = 0.0
  [../]
  [./y_top]
    type = DirichletBC
    variable = disp_y
    boundary = 'block_top'
    value = 0.0
  [../]

  #[./ns-1_bot]
  #  type = DirichletBC
  #  variable = disp_x
  #  boundary = NS-1
  #  value = 0.0
  #[../]
[]

[Materials]

  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    block = 'block_NMC'
    youngs_modulus = 177.5e3
    poissons_ratio = 0.33
  [../]
  [./isotropic_plasticity]
    type = IsotropicPlasticityStressUpdate
    block = 'block_NMC'
    yield_stress = 12600.0
    hardening_function = hf
  [../]
  [./small_stress]
    type = ComputeFiniteStrainElasticStress
    block = 'block_NMC'
  [../]
  [./thermal_expansion_strain]
    type = ComputeThermalExpansionEigenstrain
    block = 'block_NMC'
    stress_free_temperature = 0
    thermal_expansion_coeff = 0.0015
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
  [./small_stress1]
    type = ComputeFiniteStrainElasticStress
    #eigenstrain_name = eigenstrain
    block = 'block_LPS'
  [../]
  [./isotropic_plasticity1]
    type = IsotropicPlasticityStressUpdate
    block = 'block_LPS'
    yield_stress = 2000.0
    hardening_function = hf1
    #eigenstrain_name = eigenstrain
  [../]
  #[./radial_return_stress1]
  #  type = ComputeMultipleInelasticStress
  #  tangent_operator = elastic
  #  inelastic_models = 'isotropic_plasticity'
  #  #eigenstrain_name = eigenstrain
  #  block = 'block_LPS'
  #[../]
  #[./eigen]
  #type = ComputeEigenstrain
  #eigenstrain_name = eigenstrain
  #eigen_base = '1e-3 1e-3 1e-3 0 0 0'
  #block = 'block_LPS'
  #[../]

  [./thermal_expansion_strain1]
    type = ComputeThermalExpansionEigenstrain
    block = 'block_LPS'
    stress_free_temperature = 0
    thermal_expansion_coeff = 0.0
    temperature = temp
    eigenstrain_name = eigenstrain
  [../]

[]

[Executioner]
  type = Transient
  solve_type = NEWTON

  l_max_its = 50
  nl_max_its = 50
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-10
  l_tol = 1e-9

  start_time = -0.05
  n_startup_steps = 1
  end_time = 1
  dt = 0.05
  dtmin = 0.001
[]

[Outputs]
  exodus = true
[]

#[ThermalContact]
#  [thermal_contact]
#    type = GapHeatTransfer
#    variable = temp
#    primary = interface
#    secondary = 8
#    emissivity_primary = 0
#    emissivity_secondary = 0
#    gap_conductance = 1.0e9
#  []
#[]


[Postprocessors]
  [./eigenstrain_xx]
    type = ElementAverageValue
    variable = eigenstrain_xx
    execute_on = 'initial timestep_end'
  [../]
  [./eigenstrain_yy]
    type = ElementAverageValue
    variable = eigenstrain_yy
    execute_on = 'initial timestep_end'
  [../]
  [./total_strain_xx]
    type = ElementAverageValue
    variable = total_strain_xx
    execute_on = 'initial timestep_end'
  [../]
  [./total_strain_yy]
    type = ElementAverageValue
    variable = total_strain_yy
    execute_on = 'initial timestep_end'
  [../]
  [./temperature]
    type = AverageNodalVariableValue
    variable = temp
    execute_on = 'initial timestep_end'
  [../]
[]
