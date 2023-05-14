[Mesh]
 [LiGr]
   type = FileMeshGenerator
   file = data/lLi24.msh
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
   end_time = 200

[]

[Postprocessors]
  [cLi_s1]
    type = ElementAverageValue
    variable = cLi
    block = 11
    execute_on = 'initial timestep_end'
  []
  [cLi_s2]
    type = ElementAverageValue
    variable = cLi
    block = 12
    execute_on = 'initial timestep_end'
  []
  [cLi_s3]
    type = ElementAverageValue
    variable = cLi
    block = 13
    execute_on = 'initial timestep_end'
  []
  [cLi_s4]
    type = ElementAverageValue
    variable = cLi
    block = 14
    execute_on = 'initial timestep_end'
  []
  [cLi_s5]
    type = ElementAverageValue
    variable = cLi
    block = 15
    execute_on = 'initial timestep_end'
  []
  [cLi_s6]
    type = ElementAverageValue
    variable = cLi
    block = 16
    execute_on = 'initial timestep_end'
  []
#  [q_left]
#    type = ADSideDiffusiveFluxAverage
#    variable = cLi
#    boundary = Gr_btm
#    diffusivity = diffusivity
#    execute_on = 'initial timestep_end'
#  []
[]

[Outputs]
  exodus = true
  file_base = rst/lLi24
  csv = true
[]
