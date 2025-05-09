### create on 05/09/2023 by Howard Tu for the AgLi paper
## Unit system used in this code. Convert all parameters to these units !!!
# INPUT   # current density: mA/cm^2,  conductivity: mS/cm, Diffusivity: um^2/s
# INTERNAL# length: um,                potential: mV,       current: nA,         time: s

[Mesh]
# uniform_refine = 1
# coord_type = RZ
 [importMesh]
   type = FileMeshGenerator
   file = data/1Ag.msh
 []
[]

[Variables]
  [potLi]
  []
[]

[AuxVariables]
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
  [ionic_conduction]
    type = ChargedTransport
    variable = potLi
    diffusivity = ionic_conductivity
    block = "blockSE"
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
  [ag_BV]
    type = ButlerVolmerIonics
    variable = potLi
    boundary = ag
#    LiPotRef = 0
    LiPotRef = -2
    ex_current= 1.3
  []
  [li_BV]
    type = ButlerVolmerIonics
    variable = potLi
    boundary = 'liMetal'
    LiPotRef = 0
    ex_current= 1.3
  []
  [cathode_current]
    type = ADNeumannBC
    variable = potLi
    boundary = blockSE_top
    value    = 0.2   ## applied current density. in unit mA/cm^2
  []

[]

[Materials/constant]
  type = Ionics
  ionic_conductivity = 0.1
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart'
  petsc_options_value = 'asm lu 1 101'

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-13
  nl_max_its = 50
  l_tol = 1e-8
[]

[VectorPostprocessors]
  [anode_current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y'
    boundary = 'ag liMetal'
    sort_by = x   
  []

  [anode_potential]
    type = SideValueSampler
    variable = 'potLi'
    boundary = 'ag liMetal'
    sort_by = x    
  []
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
  file_base = rst/1Ag
  csv = true
[]
