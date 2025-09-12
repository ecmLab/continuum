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
  [./axial]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
#  [./disp]
#    type = DirichletBC
#    variable = disp_x
#    boundary = right
#    value = 0.1
#  [../]
  [./press]
    type = ADPressure
    variable = disp_x
    boundary = right
    factor = -1e9
  [../]
  [./bottom_x]
    type = DirichletBC
    variable = disp_x 
    boundary = bottom
    value = 0
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
    boundary = bottom
    value = 0
  [../]
[]

[Materials]
  [./elasticity]
    type = ADComputeIsotropicElasticityTensor
    poissons_ratio = 0.3
    youngs_modulus = 1e10
  [../]
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

  petsc_options_iname = -pc_hypre_type
  petsc_options_value = boomeramg

  dtmin = 0.05
  num_steps = 1
[]

[Outputs]
  exodus = true
[]
