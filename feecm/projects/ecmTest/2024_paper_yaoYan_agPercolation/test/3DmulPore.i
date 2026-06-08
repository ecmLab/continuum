[Mesh]
 [GeneratedMesh]
   type = FileMeshGenerator
   file = data/mdl1.msh
 []

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
    block = blockSE
  []
  [iLi_y]
    order = CONSTANT
    family = MONOMIAL
    block = blockSE
  []
  [iLi_z]
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
  [iEn_z]
    order = CONSTANT
    family = MONOMIAL
  []
  [LiPot]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [ionic_conduction]
    type = ParamDiffusion
    variable = potLi
    conductivity = ionic_conductivity
  []

  [electronic_conduction]
    type = ParamDiffusion
    variable = potEn
    conductivity = electronic_conductivity
  []

[]

[AuxKernels]
  [iLi_x]
    type = CurrentDensityLi
    variable = iLi_x
    component = x
    execute_on = timestep_end
    potential = potLi
  []
  [iLi_y]
    type = CurrentDensityLi
    variable = iLi_y
    component = y
    execute_on = timestep_end
    potential = potLi
  []
  [iLi_z]
    type = CurrentDensityLi
    variable = iLi_z
    component = z
    execute_on = timestep_end
    potential = potLi
  []
  [iEn_x]
    type = CurrentDensityEn
    variable = iEn_x
    component = x
    execute_on = timestep_end
    potential = potEn
  []
  [iEn_y]
    type = CurrentDensityEn
    variable = iEn_y
    component = y
    execute_on = timestep_end
    potential = potEn
  []
  [iEn_z]
    type = CurrentDensityEn
    variable = iEn_z
    component = z
    execute_on = timestep_end
    potential = potEn
  []
  [LiPotInSE]
    type = LiPotInSE
    variable = LiPot
    execute_on = timestep_end
    potLi = potLi
    potEn = potEn
  []
[]

[BCs]
  [anodeLi_BV]
    type = DepositionElectrode
    variable = potLi
    potEn    = potEn
    LiPotElectrode = LiPotAnode
    boundary = left
  []
  [cathodeLi_BV]
    type = DepositionElectrode
    variable = potLi
    potEn = potEn
    LiPotElectrode = LiPotCathode
    boundary = right
  []
  [poreLi_BV]
    type = DepositionVoid
    variable = potLi
    potRef = potEn
    boundary = pore
  []

  [anodeEn_potential]
    type = DirichletBC
    variable = potEn
    boundary = left
    value = 0
  []
  [cathodeEn_current]
    type = NeumannEn
    variable = potEn
    potLi    = potLi
    LiPotElectrode = LiPotCathode
    boundary = right
  []
  [poreEn_BV]
    type = DepositionVoid
    variable = potEn
    potRef = potLi
    boundary = pore
  []

[]

[Postprocessors]
  [anode_En]
    type = SideTotEn
    variable = potEn
    boundary = left
    diffusivity = electronic_conductivity
  []
  [cathode_En]
    type = SideTotEn
    variable = potEn
    boundary = right
    diffusivity = electronic_conductivity
  []
[]

[VectorPostprocessors]
  [anode_potential]
    type = SideValueSampler
    variable = 'potLi potEn'
    boundary = left
    sort_by = y
  []
#  [side_potential]
#    type = SideValueSampler
#    variable = 'potLi potEn'
#    boundary = 'top side'
#    sort_by = x
#  []
  [poreSE_potential]
    type = SideValueSampler
    variable = 'potLi potEn'
    boundary = pore
    sort_by = y
  []
  [cathode_potential]
    type = SideValueSampler
    variable = 'potLi potEn'
    boundary = right
    sort_by = y
  []

  [anode_current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y iLi_z iEn_x iEn_y iEn_z'
    boundary = left
    sort_by = y
  []
  [poreSE_current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y iLi_z iEn_x iEn_y iEn_z'
    boundary = pore
    sort_by = y
  []
  [cathode_current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y iLi_z iEn_x iEn_y iEn_z'
    boundary = right
    sort_by = y
  []

[]

[Materials]
## Unit system used in this code:
## length: um,     potential: V,    current: nA
## current density: mA/cm^2,  conductivity: mS/cm
## Convert all parameters to these units !!!

  [packed]
    type = PackedColumn
    inlet_current = 0.2
    exchange_current = 1.3
    reaction_rate = 0.5
    ionic_conductivity = 0.2
    electronic_conductivity = 0.0001
    metal_conductivity = 1000
    LiPotAnode = 0
    LiPotCathode = 0
  []
[]

#[Preconditioning]
#  [./SMP]
#    type = SMP
#    full = true
#    solve_type = 'NEWTON'
#  [../]
#[]

[Executioner]
  type = Steady
  solve_type = NEWTON

#  nl_rel_tol = 1e-4
#  nl_abs_tol = 1e-8
#  nl_max_its = 6
#  l_tol = 1e-4
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
  file_base = rst/mdl1
  csv = true
[]
