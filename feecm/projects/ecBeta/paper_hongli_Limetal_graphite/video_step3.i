[Mesh]
# uniform_refine = 1
 [importMesh]
   type = FileMeshGenerator
   file = rst/video/video_l45.e
   use_for_exodus_restart = true
 []
  [./createNewSidesetOne]
    type = SideSetsFromBoundingBoxGenerator
    input = importMesh
    boundary_id_old = 'Gr_top'
    boundary_id_new = 101
    bottom_left = '-0.01 9.0 0'
    top_right = '22.0 10.5 0'
    block_id = 7
  []
[]

[Variables]
  [cLi]
   initial_from_file_var = cLi
   block = 'blockGr'
  []
  [cPl]
   block = 'blockLi'
  []
[]

[Functions]
  [tPl]
     type = ParsedFunction
     value = 'if (x<=22, 0.1667, 0)'
  []
[]

[ICs]
  [cIC]
    type = FunctionIC
    variable = 'cPl'
    function = tPl
  []
[]

[Kernels]
  [cLi_time]
    type = ADTimeDerivative
    variable = cLi
  []
  [cLi_diff]
    type = ADMatDiffusion
    variable = cLi
    diffusivity = diffusivity
  []

  [cPl_time]
    type = ADTimeDerivative
    variable = cPl
  []
  [cPl_diff]
    type = ADMatDiffusion
    variable = cPl
    diffusivity = diffusivity1
  []
[]

[BCs]
  [cLi_left]
    type = ADDirichletBC
    variable = cLi
    boundary = 'Li_lft Li_top'
    value = 0.1667
  []
  [cLi_top]
    type = ADDirichletBC
    variable = cLi
    boundary = 101
    value = 0.1667
  []
    [Gr_anode_BV]
    type = ADNeumannBC
    variable = cLi
    boundary = Gr_rgt
    value = 0.001
  []
  [cPl_left]
    type = ADDirichletBC
    variable = cPl
    boundary = Pl_lft
    value = 0.1667
  []
[]

[Materials]
## Unit system used in this code:
## length: um,     potential: V,    current: nA
## Diffusivity: um^2/s
## current density: mA/cm^2,  conductivity: mS/cm
## Convert all parameters to these units !!!

  [cPl]
  type = ADGenericConstantMaterial
  prop_names = 'diffusivity1'
  prop_values = 0.00001
  []
  [cLi]
  type = paperHongli
  cLi  = cLi
  c1   = 0.5
  c2   = -2
  block = 'blockGr'
 []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  dtmin = 0.01
  dt   = 10
  start_time = 2020
  end_time = 2200

[]

[Outputs]
  exodus = true
  file_base = rst/video/video_l46
  csv = true
[]
