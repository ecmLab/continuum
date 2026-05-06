name = 'L2'

[Mesh]
 [GeneratedMesh]
   type = FileMeshGenerator
   file = data/${name}.msh
 []

 [./interface]
   type = SideSetsBetweenSubdomainsGenerator
   input = GeneratedMesh
   primary_block = 'blockSE'
   paired_block = 'blockAnd'
   new_boundary = 'interface0'
 [../]

 [./interface1]
   type = SideSetsBetweenSubdomainsGenerator
   input = interface
   primary_block = 'blockSE'
   paired_block = 'blockCtd'
   new_boundary = 'interface1'
 [../]

 [./break_boundary]
   input = interface1
   type = BreakBoundaryOnSubdomainGenerator
 [../]

 [./interface2]
   type = SideSetsBetweenSubdomainsGenerator
   input = break_boundary
   primary_block = 'blockAnd'
   paired_block = 'blockSE'
   new_boundary = 'interface2'
 [../]

 [./interface3]
   type = SideSetsBetweenSubdomainsGenerator
   input = interface2
   primary_block = 'blockCtd'
   paired_block = 'blockSE'
   new_boundary = 'interface3'
 [../]

[]

[Variables]
  [potLi]
   block = 'blockSE'
  []
  [potEn]
   block = 'blockAnd blockCtd'
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
    type = ParamDiffusion
    variable = potLi
    conductivity = ionic_conductivity
    block = blockSE
  []

  [metal_conduction]
    type = ParamDiffusion
    variable = potEn
    conductivity = metal_conductivity
    block = 'blockAnd blockCtd'
  []

[]

[InterfaceKernels]
  [interface]
    type = SingleSEInterface
    variable = potLi
    neighbor_var  = potEn
    boundary = 'interface0 interface1'
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
    block = blockSE
  []
  [iLi_y]
    type = CurrentDensity
    variable = iLi_y
    component = y
    conductivity = ionic_conductivity
    execute_on = timestep_end
    potential = potLi
    block = blockSE
  []
  [iEn_x]
    type = CurrentDensity
    variable = iEn_x
    component = x
    conductivity = metal_conductivity
    execute_on = timestep_end
    potential = potEn
    block = 'blockAnd blockCtd'
  []
  [iEn_y]
    type = CurrentDensity
    variable = iEn_y
    component = y
    conductivity = metal_conductivity
    execute_on = timestep_end
    potential = potEn
    block = 'blockAnd blockCtd'
  []
[]

[BCs]
#  [SE_anode_BV]
#    type = SingleSEElectrodeBV
#    variable = potLi
#    LiPotElectrode = LiPotAnode
#    boundary = SE_left
#  []
#  [SE_cathode_current]
#    type = SingleSEElectrodeNeumann
#    variable = potLi
#    boundary = SE_right
#  []

  [Li_anode_potential]
    type = DirichletBC
    variable = potEn
    boundary = And_left
    value = 0
  []

  [cathode_current]
    type = SingleSEElectrodeNeumann
    variable = potEn
    boundary = Ctd_right
  []

[]

[Postprocessors]
  [AndLft_EnTotCurrent]
    type = SideTotCrnt
    variable = potEn
    boundary = And_left
    conductivity = metal_conductivity
  []
  [SEAndInt_LiTotCurrent]
    type = SideTotCrnt
    variable = potLi
    boundary = interface0
    conductivity = ionic_conductivity
  []
  [SECtdInt_LiTotCurrent]
    type = SideTotCrnt
    variable = potLi
    boundary = interface1
    conductivity = ionic_conductivity
  []
  [CtdRgt_EnTotCurrent]
    type = SideTotCrnt
    variable = potEn
    boundary = Ctd_right
    conductivity = metal_conductivity
  []
[]

[VectorPostprocessors]
  [AndTop_potEn]
    type = SideValueSampler
    variable = 'potEn'
    boundary = And_top
    sort_by = x
  []
  [SEAndInt_potEn]
    type = SideValueSampler
    variable = 'potEn'
    boundary = interface2
    sort_by = x
  []
  [SEAndInt_potLi]
    type = SideValueSampler
    variable = 'potLi'
    boundary = interface0
    sort_by = x
  []
  [SETop_potLi]
    type = SideValueSampler
    variable = 'potLi'
    boundary = SE_top
    sort_by = x
  []
  [SEBtm_potLi]
    type = SideValueSampler
    variable = 'potLi'
    boundary = SE_bottom
    sort_by = x
  []
  [SECtdInt_potLi]
    type = SideValueSampler
    variable = 'potLi'
    boundary = interface1
    sort_by = y
  []
  [SECtdInt_potEn]
    type = SideValueSampler
    variable = 'potEn'
    boundary = interface3
    sort_by = y
  []
  [CtdRgt_potEn]
    type = SideValueSampler
    variable = 'potEn'
    boundary = Ctd_right
    sort_by = y
  []

  [AndLft_crntEn]
    type = SideValueSampler
    variable = 'iEn_x iEn_y'
    boundary = And_left
    sort_by = y
  []
  [SEAndInt_crntEn]
    type = SideValueSampler
    variable = 'iEn_x iEn_y'
    boundary = interface2
    sort_by = x
  []
  [SEAndInt_crntLi]
    type = SideValueSampler
    variable = 'iLi_x iLi_y'
    boundary = interface0
    sort_by = x
  []
  [SECtdInt_crntLi]
    type = SideValueSampler
    variable = 'iLi_x iLi_y'
    boundary = interface1
    sort_by = y
  []
  [SECtdInt_crntEn]
    type = SideValueSampler
    variable = 'iEn_x iEn_y'
    boundary = interface3
    sort_by = y
  []
  [CtdRgt_crntEn]
    type = SideValueSampler
    variable = 'iEn_x iEn_y'
    boundary = Ctd_right
    sort_by = y
  []

[]

[Materials]
## Unit system used in this code:
## length: um,     potential: V,    current: nA
## current density: mA/cm^2,  conductivity: mS/cm
## if applied_current < 0, Li+ flowes into the SE from the boudanry
## Convert all parameters to these units !!!

  [packed]
    type = SingleSE
    applied_current = -1.0
    exchange_current = 1.3
    reaction_rate = 0.5
    ionic_conductivity = 1
    gb_conductivity = 0
    metal_conductivity = 100000
    LiPotAnode = 0
    LiPotCathode = 0
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart'
  petsc_options_value = 'asm lu 1 101'

#  line_search = 'none'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-10
  l_tol = 1e-6
  nl_max_its = 50

[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
  file_base = rst/${name}
  csv = true
[]
