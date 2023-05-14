[Mesh]
 uniform_refine = 1
 [importMesh]
   type = FileMeshGenerator
   file = rst/lLi24.e
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
    LiPotElectrode = 0.0
    boundary = Gr_Li
  []
  [Gr_SE_current]
    type = SingleSEElectrodeNeumann
    variable = potLi
    boundary = Gr_rgt
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
  applied_current = 0.2
  exchange_current = 1.3
  reaction_rate = 0.5
  ionic_conductivity = 0.1
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

  dtmin = 0.01
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    growth_factor = 1.5
    cutback_factor = 0.25
    optimal_iterations = 10
  [../]
   end_time = 10

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
    boundary = Gr_btm
    diffusivity = diffusivity
    execute_on = 'initial timestep_end'
  []
[]

[Outputs]
  exodus = true
  file_base = rst/lLi24_ec
  csv = true
[]
