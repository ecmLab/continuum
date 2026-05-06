[Mesh]
 [LiGr]
   type = FileMeshGenerator
   file = data/t24.msh
 []

[section1]
   type = ParsedSubdomainMeshGenerator
   input = LiGr
   combinatorial_geometry = 'x <= 4.425'
   block_id = 11
 []
 [section2]
   type = ParsedSubdomainMeshGenerator
   input = section1
   combinatorial_geometry = 'x > 4.425 & x <= 8.85'
   block_id = 12
 []
 [section3]
   type = ParsedSubdomainMeshGenerator
   input = section2
   combinatorial_geometry = 'x > 8.85 & x <= 13.275'
   block_id = 13
 []
 [section4]
   type = ParsedSubdomainMeshGenerator
   input = section3
   combinatorial_geometry = 'x > 13.275 & x <= 17.7'
   block_id = 14
 []
 [section5]
   type = ParsedSubdomainMeshGenerator
   input = section4
   combinatorial_geometry = 'x > 17.7 & x <= 22.125'
   block_id = 15
 []
 [section6]
   type = ParsedSubdomainMeshGenerator
   input = section5
   combinatorial_geometry = 'x > 22.125'
   block_id = 16
 []

[]

[Variables]
  [./cLi]
    order = THIRD
    family = HERMITE
    initial_condition = 0.1
  [../]
[]

[Kernels]
  [./ie_c]
    type = ADTimeDerivative
    variable = cLi
  [../]
  [./CHSolid]
    type = CahnHilliard
    variable = cLi
    f_name = F
    mob_name = M
  [../]
#  [./CHInterface]
#    type = CHInterface
#    variable = cLi
#    mob_name = M
#    kappa_name = kappa_c
#  [../]
[]

#[BCs]
#  [left]
#    type = ADDirichletBC
#    variable = cLi
#    boundary = Gr_Li
#    value = 0.1667
#  []
#[]

[Materials]
  [./consts]
    type = GenericConstantMaterial
    prop_names  = 'M kappa_c'
    prop_values = '0.00000001 0.01'
  [../]
  [./free_energy]
    type = DerivativeParsedMaterial
    f_name = F
    args = 'cLi'
    function = '(1-cLi)^2 * (1+cLi)^2'
  [../]
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 31'
  dtmin = 0.01
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    growth_factor = 1.5
    cutback_factor = 0.25
    optimal_iterations = 10
  [../]
   end_time = 200

[]

[Outputs]
  exodus = true
  file_base = rst/t24
  csv = true
[]
