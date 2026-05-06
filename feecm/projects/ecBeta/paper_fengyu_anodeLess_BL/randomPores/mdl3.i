[Mesh]
 [GeneratedMesh]
   type = FileMeshGenerator
   file = data/mdl3.msh
 []

 [./interface]
   type = SideSetsBetweenSubdomainsGenerator
   input = GeneratedMesh
   master_block = 'blockSE'
   paired_block = 'blockLi'
   new_boundary = 'interface01'
 [../]

 [./break_boundary]
   input = interface
   type = BreakBoundaryOnSubdomainGenerator
 [../]

[]

[Variables]
  [potLi]
   block = 'blockSE'
  []
  [potEn]
   block = 'blockSE'
  []
  [potMt]
   block = 'blockLi'
  []
[]

[AuxVariables]
  [iLi_x]
    order = CONSTANT
    family = MONOMIAL
    block = blockSE
  []
  [iLi_y]
    order = CONSTANT
    family = MONOMIAL
    block = blockSE
  []
  [iEn_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [iEn_y]
    order = CONSTANT
    family = MONOMIAL
  []
 [iMt_x]
    order = CONSTANT
    family = MONOMIAL
    block = blockLi
  []
  [iMt_y]
    order = CONSTANT
    family = MONOMIAL
    block = blockLi
  []
  [LiPot]
    order = CONSTANT
    family = MONOMIAL
    block = blockSE
  []
[]

[Kernels]
  [ionic_conduction]
    type = ParamDiffusion
    variable = potLi
    conductivity = ionic_conductivity
    block = blockSE
  []

  [electronic_conduction]
    type = ParamDiffusion
    variable = potEn
    conductivity = electronic_conductivity
    block = blockSE
  []

  [metal_conduction]
    type = ParamDiffusion
    variable = potMt
    conductivity = metal_conductivity
    block = blockLi
  []
[]

[InterfaceKernels]
  [interface]
    type = MixSEInterfaceEn
    variable = potEn
    neighbor_var = potMt
    potLi    = potLi
    boundary = interface01
  []
  [interface_potential]
    type = PenaltyInterfaceDiffusion
    variable = potEn
    neighbor_var = potMt
    boundary = interface01
    penalty = 1e6
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
    conductivity = electronic_conductivity
    execute_on = timestep_end
    potential = potEn
    block = blockSE
  []
  [iEn_y]
    type = CurrentDensity
    variable = iEn_y
    component = y
    conductivity = electronic_conductivity
    execute_on = timestep_end
    potential = potEn
    block = blockSE
  []
  [iMt_x]
    type = CurrentDensity
    variable = iMt_x
    component = x
    conductivity = metal_conductivity
    execute_on = timestep_end
    potential = potMt
    block = blockLi
  []
  [iMt_y]
    type = CurrentDensity
    variable = iMt_y
    component = y
    conductivity = metal_conductivity
    execute_on = timestep_end
    potential = potMt
    block = blockLi
  []
  [LiPotInSE]
    type = LiPotInSE
    variable = LiPot
    execute_on = timestep_end
    potLi = potLi
    potEn = potEn
    block = blockSE
  []
[]

[BCs]
  [anodeLi_BV]
    type = MixSEIntBV
    variable = potLi
    potEn    = potEn
    LiPotElectrode = LiPotAnode
    boundary = top
  []
  [cathodeLi_BV]
    type = MixSEIntBV
    variable = potLi
    potEn = potEn
    LiPotElectrode = LiPotCathode
    boundary = bottom
  []
  [poreLi_BV]
    type = MixSEVoidBV
    variable = potLi
    potRef = potEn
    boundary = pore
  []

  [anodeEn_potential]
    type = DirichletBC
    variable = potEn
    boundary = top
    value = 0
  []
  [cathodeEn_current]
    type = NeumannEn
    variable = potEn
    potLi    = potLi
    LiPotElectrode = LiPotCathode
    boundary = bottom
  []
  [poreEn_BV]
    type = MixSEVoidBV
    variable = potEn
    potRef = potLi
    boundary = pore
 []

[]

[Postprocessors]
  [anode_En]
    type = SideTotCrnt
    variable = potEn
    boundary = top
    conductivity = electronic_conductivity
  []
  [cathode_En]
    type = SideTotCrnt
    variable = potEn
    boundary = bottom
    conductivity = electronic_conductivity
  []
[]

[VectorPostprocessors]
  [anode_potential]
    type = SideValueSampler
    variable = 'potLi potEn'
    boundary = top
    sort_by = y
  []
  [cathode_potential]
    type = SideValueSampler
    variable = 'potLi potEn'
    boundary = bottom
    sort_by = y
  []

  [anode_current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y iEn_x iEn_y'
    boundary = top
    sort_by = y
  []
  [cathode_current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y iEn_x iEn_y'
    boundary = bottom
    sort_by = y
  []

[]

[Materials]
## Unit system used in this code:
## length: um,     potential: V,    current: nA
## current density: mA/cm^2,  conductivity: mS/cm
## Convert all parameters to these units !!!

 [packed]
    type = MixSE
    inlet_current = 0.2
    exchange_current = 1.3
    reaction_rate = 0.5
    ionic_conductivity = 0.3
    metal_conductivity = 10
    electronic_conductivity = 0.0001
    LiPotAnode = 0
    LiPotCathode = 0
    electron_concentration = 0.008
  []

[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart'
  petsc_options_value = 'asm lu 1 101'

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-10
  nl_max_its = 50
  l_tol = 1e-6
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
  file_base = rst/mdl3
  csv = true
[]
