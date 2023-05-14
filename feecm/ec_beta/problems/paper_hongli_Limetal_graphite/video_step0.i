[Mesh]
 [LiGr]
   type = FileMeshGenerator
   file = data/video.msh
 []
[]

[Variables]
  [cLi]
   initial_condition = 0
   block = 'blockGr'
  []
  [cPl]
   block = 'blockLi'
  []
[]

[Functions]
  [tPl]
     type = ParsedFunction
     value = 'if (x-0.1*t-0.5<0, 0.1667, 0)'
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
  [ctLi_left]
    type = ADDirichletBC
    variable = cLi
    boundary = 'Li_lft Li_top'
    value = 0.1667
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
  prop_values = 0.0001
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
  dt   = 2.5
  end_time = 5

[]

[Outputs]
  exodus = true
  file_base = rst/video/video_l1
  csv = true
[]
