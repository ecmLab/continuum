## create on 05/09/2023 by Howard Tu for the AgLi paper
## Unit system used in this code. Convert all parameters to these units !!!
# INPUT   # current density: mA/cm^2,  conductivity: mS/cm, Diffusivity: um^2/s
# INTERNAL# length: um,                potential: mV,       current: nA,         time: s

[Mesh]
# uniform_refine = 1
 [importMesh]
   type = FileMeshGenerator
   file = data/purePore.msh
 []
[]

[Variables]
  [potLi]
  []
  [potEn]
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
  [iEn_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [iEn_y]
    order = CONSTANT
    family = MONOMIAL
  []

[]

[Kernels]
  [ionic_conduction]
    type = ChargedTransport
    variable = potLi
    diffusivity = ionic_conductivity
    block = "blockBL"
  []

  [electronic_conduction]
    type = ChargedTransport
    variable = potEn
    diffusivity = electronic_conductivity
    block = "blockBL"
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

  [iEn_x]
    type = CurrentDensity
    variable = iEn_x
    component = x
    conductivity = electronic_conductivity
    execute_on = timestep_end
    potential = potEn
  []
  [iEn_y]
    type = CurrentDensity
    variable = iEn_y
    component = y
    conductivity = electronic_conductivity
    execute_on = timestep_end
    potential = potEn
  []
[]

[BCs]
  [delithiation_BV]    ##Lithiation at the SE/BL interface
    type = ButlerVolmerMiec
    variable = potLi
    potEn = potEn
    boundary = blockBL_top
    LiCrtRef = -0.68
    LiPotRef = 0.0
    ex_current= 1.3
  []
  [Ag_BV]            ## Li-Ag alloying at the Ag/C interface inside the BL
    type = ButlerVolmerMiec
    variable = potLi
    potEn = potEn
    boundary = pore
    LiCrtRef = 0.0
    LiPotRef = -10
    ex_current= 1.3
  []
  [lithiation_BV]      ## Lithiation at the CC/BL interface
    type = ButlerVolmerMiec
    variable = potLi
    potEn = potEn
    boundary = blockBL_btm
    LiCrtRef = 0.0
    LiPotRef = 0.0
    ex_current= 1.3
  []

  [coupling_current_SE]   ## The coupling of electronic and ionic current at the SE/BL interface
    type = CouplingMiec
    variable = potEn
    potLi  = potLi
    boundary = blockBL_top
    LiCrtRef  = -0.68   ## reference current at this interface, in unit mA/cm^2; negative means cell charging
  []
  [coupling_current_pore]   ## The coupling of electronic and ionic current at the pore surface
    type = CouplingMiec
    variable = potEn
    potLi  = potLi
    boundary = pore
    LiCrtRef  = 0.0   ## reference current at this interface, in unit mA/cm^2
  []
  [coupling_current_CC]   ## The coupling of electronic and ionic current at the CC/BL interface
    type = CouplingMiec
    variable = potEn
    potLi  = potLi
    boundary = blockBL_btm
    LiCrtRef  = 0.68   ## reference current at this interace, in unit mA/cm^2; 
  []

  [anodeEn_potential]
    type = DirichletBC
    variable = potEn
    boundary = blockBL_btm
    value = 0
  []

[]

[Materials/constant]
  type = Miec
  ionic_conductivity = 10
  electronic_conductivity = 1000
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

#[VectorPostprocessors]
#  [anode_current]
#    type = SideValueSampler
#    variable = 'iLi_x iLi_y'
#    boundary = 'blockSE_btm_lft interface blockSE_btm_rgt'
#    sort_by = x
#  []

#  [anode_potential]
#    type = SideValueSampler
#    variable = 'potLi'
#    boundary = 'blockSE_btm_lft interface blockSE_btm_rgt'
#    sort_by = x
#  []
#[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
  file_base = rst/mdl
  csv = true
[]
