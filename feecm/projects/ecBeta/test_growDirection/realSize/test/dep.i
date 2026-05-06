name = 'mdl2'

[Mesh]
 [GeneratedMesh]
   type = FileMeshGenerator
   file = data/${name}.msh
 []

 [./interface]
   type = SideSetsBetweenSubdomainsGenerator
   input = GeneratedMesh
   primary_block = 'blockSE'
   paired_block = 'blockLi'
   new_boundary = 'interface0'
 [../]

 [./break_boundary]
   input = interface
   type = BreakBoundaryOnSubdomainGenerator
 [../]

 [./interface1]
   type = SideSetsBetweenSubdomainsGenerator
   input = break_boundary
   primary_block = 'blockLi'
   paired_block = 'blockSE'
   new_boundary = 'interface1'
 [../]

[]

[Variables]
  [potLi]
   block = 'blockSE'
  []
  [potEn]
   block = 'blockLi'
  []
[]

[AuxVariables]
  [iLi_x]
    order = CONSTANT
    family = MONOMIAL
#    block = blockSE
  []
  [iLi_y]
    order = CONSTANT
    family = MONOMIAL
#    block = blockSE
  []
  [iEn_x]
    order = CONSTANT
    family = MONOMIAL
#    block = blockLi
  []
  [iEn_y]
    order = CONSTANT
    family = MONOMIAL
#    block = blockLi
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
    block = blockLi
  []

[]

[InterfaceKernels]
  [interface]
    type = SingleSEInterface
    variable = potLi
    neighbor_var  = potEn
    boundary = interface0
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
    block = blockLi
  []
  [iEn_y]
    type = CurrentDensity
    variable = iEn_y
    component = y
    conductivity = metal_conductivity
    execute_on = timestep_end
    potential = potEn
    block = blockLi
  []
[]

[BCs]
#  [SE_anode_BV]
#    type = SingleSEElectrodeBV
#    variable = potLi
#    LiPotElectrode = LiPotAnode
#    boundary = SE_left
#  []
  [SE_cathode_current]
    type = SingleSEElectrodeNeumann
    variable = potLi
    boundary = SE_right
  []

  [Li_anode_potential]
    type = DirichletBC
    variable = potEn
    boundary = Li_left
    value = 0
  []

[]

[Postprocessors]
#  [SEAnode_TotCurrent]
#    type = SideTotCrnt
#    variable = potLi
#    boundary = SE_left
#    conductivity = ionic_conductivity
#  []
  [SEInt_TotCurrent]
    type = SideTotCrnt
    variable = potLi
    boundary = interface0
    conductivity = ionic_conductivity
  []
  [SECathode_TotCurrent]
    type = SideTotCrnt
    variable = potLi
    boundary = SE_right
    conductivity = ionic_conductivity
  []
  [LiAnode_TotCurrent]
    type = SideTotCrnt
    variable = potEn
    boundary = Li_left
    conductivity = metal_conductivity
  []
[]

[VectorPostprocessors]
#  [SEAnode_potential]
#    type = SideValueSampler
#    variable = 'potLi'
#    boundary = SE_left
#    sort_by = y
#  []
  [SEInt_potential]
    type = SideValueSampler
    variable = 'potLi'
    boundary = interface0
    sort_by = x
  []
  [SETop_potential]
    type = SideValueSampler
    variable = 'potLi'
    boundary = SE_top
    sort_by = x
  []
  [SECathode_potential]
    type = SideValueSampler
    variable = 'potLi'
    boundary = SE_right
    sort_by = y
  []
  [SEBottom_potential]
    type = SideValueSampler
    variable = 'potLi'
    boundary = SE_bottom
    sort_by = x
  []
#  [SEAnode_current]
#    type = SideValueSampler
#    variable = 'iLi_x iLi_y'
#    boundary = SE_left
#    sort_by = y
#  []
  [SEInt_current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y'
    boundary = interface0
    sort_by = x
  []
  [SECathode_current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y'
    boundary = SE_right
    sort_by = y
  []

  [LiInt_potential]
    type = SideValueSampler
    variable = 'potEn'
    boundary = interface1
    sort_by = x
  []
  [LiCathode_potential]
    type = SideValueSampler
    variable = 'potEn'
    boundary = Li_right
    sort_by = y
  []
  [LiTop_potential]
    type = SideValueSampler
    variable = 'potEn'
    boundary = Li_top
    sort_by = x
  []
  [LiAnode_current]
    type = SideValueSampler
    variable = 'iEn_x iEn_y'
    boundary = Li_left
    sort_by = x
  []
  [LiInt_current]
    type = SideValueSampler
    variable = 'iEn_x iEn_y'
    boundary = interface1
    sort_by = x
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

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-10
  l_tol = 1e-6
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
  file_base = rst/${name}
  csv = true
[]
