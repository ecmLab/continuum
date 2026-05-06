[Mesh]
# uniform_refine = 1
 [importMesh]
   type = FileMeshGenerator
   file = rst/video/vid_l45.e
   use_for_exodus_restart = true
 []

[section1]
   type = ParsedSubdomainMeshGenerator
   input = importMesh 
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
   initial_from_file_var = cLi
  []
[]

[AuxVariables]
  [./bounds_dummy]
    order = FIRST
    family = LAGRANGE
  [../]
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
  [Gr_anode_BV]
    type = ADNeumannBC
    variable = cLi
    boundary = Gr_rgt
    value = 0.01
  []
[]

[Bounds]
  [./cLi_bound]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = cLi
    bound_type = upper
    bound_value = 0.1667
  [../]
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
#  solve_type = 'PJFNK'
  petsc_options_iname = '-snes_max_linear_solve_fail -ksp_max_it -pc_type -sub_pc_factor_levels -snes_linesearch_type -snes_type'
  petsc_options_value = '0                           30          asm      16                    basic                 vinewtonrsls'

  dtmin = 0.1
#  [./TimeStepper]
#    type = IterationAdaptiveDT
#    dt = 10
#    growth_factor = 1.5
#    cutback_factor = 0.25
#    optimal_iterations = 10
#  [../]
  dt = 2.5
  start_time = 2000
  end_time   = 2180

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
[]

[Outputs]
  exodus = true
  file_base = rst/video/vid_l46
  csv = true
[]
