[Mesh]
# uniform_refine = 1
 [importMesh]
   type = FileMeshGenerator
   file = rst/t24.e
   use_for_exodus_restart = true
 []
[]

[Variables]
  [cLi]
   initial_from_file_var = cLi
  []
  [potLi]
   initial_condition = 0
  []
[]

[AuxVariables]
  [./bounds_dummy]
    order = FIRST
    family = LAGRANGE
  [../]

  [iLi_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [iLi_y]
    order = CONSTANT
    family = MONOMIAL
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

  [ionic_conduction]
    type = ParamDiffusion
    variable = potLi
    conductivity = ionic_conductivity
  []

[]

[AuxKernels]
  [iLi_x]
    type = CurrentDensity
    variable = iLi_x
    component = x
    conductivity = ionic_conductivity
    execute_on = timestep_end
    potential = potLi
  []
  [iLi_y]
    type = CurrentDensity
    variable = iLi_y
    component = y
    conductivity = ionic_conductivity
    execute_on = timestep_end
    potential = potLi
  []
[]

[BCs]
  [Li_anode_potential]
    type = ADDirichletBC
    variable = potLi
    boundary = Gr_Li
    value = 0.0
  []
  [Gr_anode_BV]
    type = SingleSEElectrodeBV
    variable = potLi
    boundary = Gr_Li
    exchange_current = 13
  []
  [Gr_SE_current]
    type = ADNeumannBC
    variable = potLi
    boundary = Gr_rgt
    value = 0.1
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
#  type = Steady
  type = Transient
  solve_type = NEWTON
#  solve_type = 'PJFNK'
  petsc_options_iname = '-snes_max_linear_solve_fail -ksp_max_it -pc_type -sub_pc_factor_levels -snes_linesearch_type -snes_type'
  petsc_options_value = '0                           30          asm      16                    basic                 vinewtonrsls'

  dtmin = 0.1
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    growth_factor = 1.5
    cutback_factor = 0.25
    optimal_iterations = 10
  [../]
   end_time = 1

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

[VectorPostprocessors]
  [current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y'
    boundary = 'Gr_Li'
    sort_by = x
    execute_on = 'timestep_end'
  []
  [potLi]
    type = NodalValueSampler
    variable = 'potLi'
    sort_by = x
    execute_on = 'timestep_end'
  []
#  [Gr_potLi]
#    type = LineValueSampler
#    start_point = '0  0'
#    eend_point  = '26.55 0'
#    variable = 'potLi'
#    num_points = 100
#    sort_by = x
#  []

[]

[Outputs]
  exodus = true
  file_base = rst/t24_miec
  csv = true
[]
