[Mesh]
 [GeneratedMesh]
   type = FileMeshGenerator
   file = data/paper_rgtSide.msh
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
    boundary = Gr_Li
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
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    growth_factor = 1.5
    cutback_factor = 0.25
    optimal_iterations = 10
  [../]
   end_time = 1800

[]

[Postprocessors]
  [cLi_avg]
    type = ElementAverageValue
    variable = cLi
    execute_on = 'initial timestep_end'
  []
  [q_left]
    type = ADSideDiffusiveFluxAverage
    variable = cLi
    boundary = Gr_SE
    diffusivity = diffusivity
    execute_on = 'initial timestep_end'
  []
[]

[Outputs]
  exodus = true
  file_base = rst/paper_rgtSide
  csv = true
[]
