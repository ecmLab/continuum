[Mesh]
# uniform_refine = 1
 [importMesh]
   type = FileMeshGenerator
   file = rst/video/video_input.e
   use_for_exodus_restart = true
 []

  [./createNewSidesetOne]
    type = SideSetsFromBoundingBoxGenerator
    input = importMesh
    boundary_id_old = 'Gr_top'
    boundary_id_new = 101
    bottom_left = '-0.01 9.0 0'
    top_right = 'liLength 10.5 0'
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
  [potLi]
   initial_condition = 0
#   block = 'blockGr'
  []
[]

[Functions]
  [tPl]
     type = ParsedFunction
     value = 'if (x-0.005*t-11<0, 0.1667, 0)'
  []
[]

[ICs]
  [cIC]
    type = FunctionIC
    variable = 'cPl'
    function = tPl
  []
[]

[AuxVariables]
  [./bounds_dummy]
    order = FIRST
    family = LAGRANGE
  [../]

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

  [ionic_conduction]
    type = ParamDiffusion
    variable = potLi
    conductivity = ionic_conductivity
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
  [cPl_left]
    type = ADDirichletBC
    variable = cPl
    boundary = Pl_lft
    value = 0.1667
  []

  [Li_lft_potential]
    type = ADDirichletBC
    variable = potLi
    boundary = 'Li_lft Li_top'
    value = 0.0
  []
  [Li_top_potential]
    type = ADDirichletBC
    variable = potLi
    boundary = 101
    value = 0.0
  []
  [Li_lft_BV]
    type = SingleSEElectrodeBV
    variable = potLi
    boundary = 'Li_lft Li_top'
    exchange_current = 13
    block = 'blockGr'
  []
  [Li_top_BV]
    type = SingleSEElectrodeBV
    variable = potLi
    boundary = 101
    exchange_current = 13
    block = 'blockGr'
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
  block = 'blockLi'
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
#  solve_type = 'PJFNK'
  petsc_options_iname = '-snes_max_linear_solve_fail -ksp_max_it -pc_type -sub_pc_factor_levels -snes_linesearch_type -snes_type'
  petsc_options_value = '0                           30          asm      16                    basic                 vinewtonrsls'
 
  dtmin = 0.01
  dt   = 25
  start_time = 0
  end_time = 5
[]

[Outputs]
  exodus = true
  file_base = rst/video/video_output
  csv = true
[]
