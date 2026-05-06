[Mesh]
 [LiGr]
   type = FileMeshGenerator
   file = data/vid.msh
 []
[]

[Variables]
  [cLi]
   initial_condition = 0
  []
[]

[Kernels]
  [time]
    type = ADTimeDerivative
    variable = cLi
  []
  [diff]
    type = ADMatDiffusion
    variable = cLi
    diffusivity = diffusivity
  []
[]

[BCs]
  [left]
    type = ADDirichletBC
    variable = cLi
    boundary = 'Li_lft Li_top'
    value = 0.1667
  []
[]

[Materials/constant]
## Unit system used in this code:
## length: um,     potential: V,    current: nA                                                                    
## Diffusivity: um^2/s
## current density: mA/cm^2,  conductivity: mS/cm
## Convert all parameters to these units !!!

#  type = ADGenericConstantMaterial
#  prop_names = 'diffusivity'
#  prop_values = 0.01

  type = paperHongli
  cLi  = cLi
  c1   = 0.5
  c2   = -2

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
  file_base = rst/video/vid_l1
  csv = true
[]
