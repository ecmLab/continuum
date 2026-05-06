[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 30
  ny = 30
  nz = 30
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  # scale with one over Young's modulus
  [./disp_x]
    scaling = 1e-10
  [../]
  [./disp_y]
    scaling = 1e-10
  [../]
  [./disp_z]
    scaling = 1e-10
  [../]
  [T]
    initial_condition = 300.0
  []
[]

[AuxVariables]
 [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
 [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
 [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
#Mechanics
  [./stress_x]
    type = ADStressDivergenceTensors
    component = 0
    variable = disp_x
  [../]
  [./stress_y]
    type = ADStressDivergenceTensors
    component = 1
    variable = disp_y
  [../]
  [./stress_z]
    type = ADStressDivergenceTensors
    component = 2
    variable = disp_z
  [../]
#Thermo
 [heat_conduction]
    type = HeatConduction
    variable = T
  []
  [time_derivative]
    type = HeatConductionTimeDerivative
    variable = T
  []
  [heat_source]
    type = HeatSource
    variable = T
    value = 5e4
  []
[]

[AuxKernels]
 [./stress_xx]
    type = ADRankTwoAux
    variable = stress_xx
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
  [../]
 [./stress_yy]
    type = ADRankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
  [../]
 [./stress_zz]
    type = ADRankTwoAux
    variable = stress_zz
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
  [../]
[]

[BCs]
#Mechanics
  [./axial]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
  [./disp]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = '0.02*t'
  [../]
  [./bottom_y]
    type = DirichletBC
    variable = disp_y 
    boundary = bottom
    value = 0
  [../]
  [./bottom_z]
    type = DirichletBC
    variable = disp_z 
    boundary = back
    value = 0
  [../]
#Thermo
  [t_left]
    type = DirichletBC
    variable = T
    value = 300
    boundary = left
  []
  [t_right]
    type = FunctionDirichletBC
    variable = T
    function = '300+5*t'
    boundary = right
  []
[]

[Materials]
  [./elasticity]
    type = ADComputeIsotropicElasticityTensor
    poissons_ratio = 0.3
    youngs_modulus = 1e10
  [../]
[thermal]
    type = HeatConductionMaterial
    thermal_conductivity = 45.0
    specific_heat = 0.5
  []
  [density]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = 8000.0
  []
 [expansion1]
    type = ComputeThermalExpansionEigenstrain
    temperature = T
    thermal_expansion_coeff = 0.001
    stress_free_temperature = 300
    eigenstrain_name = thermal_expansion
  []
[]

[Materials]
  [./strain]
    type = ADComputeSmallStrain
  [../]
  [./stress]
    type = ADComputeLinearElasticStress
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
  dt = 0.05
  solve_type = 'NEWTON'

  petsc_options_iname = -pc_type
  petsc_options_value = lu

  dtmin = 0.05
  end_time = 5
[]

[Outputs]
  exodus = true
[]
