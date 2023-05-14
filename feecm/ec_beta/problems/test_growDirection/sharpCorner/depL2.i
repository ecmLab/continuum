[Mesh]
 [GeneratedMesh]
   type = FileMeshGenerator
   file = data/mdlL2.msh
 []

 [./interface]
   type = SideSetsBetweenSubdomainsGenerator
   input = GeneratedMesh
   master_block = 'blockSE'
   paired_block = 'blockLi'
   new_boundary = 'interface0'
 [../]

 [./break_boundary]
   input = interface
   type = BreakBoundaryOnSubdomainGenerator
 [../]

[]

[NodalNormals]
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
    block = blockLi
  []

[]

[InterfaceKernels]
  [interface]
    type = SingleSEInterfaceBV
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
  [SE_anode_BV]
    type = SingleSEElectrodeBV
    variable = potLi
    LiPotElectrode = LiPotAnode
    boundary = SE_left
  []
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
  [SEAnode_TotCurrent]
    type = SideTotLi
    variable = potLi
    boundary = SE_left
    diffusivity = ionic_conductivity
  []
  [SEInt_TotCurrent]
    type = SideTotLi
    variable = potLi
    boundary = interface0
    diffusivity = ionic_conductivity
  []
  [SECathode_TotCurrent]
    type = SideTotLi
    variable = potLi
    boundary = SE_right
    diffusivity = ionic_conductivity
  []
  [LiAnode_TotCurrent]
    type = SideTotEn
    variable = potEn
    boundary = Li_left
    diffusivity = metal_conductivity
  []
[]

[VectorPostprocessors]
  [SEAnode_potential]
    type = SideValueSampler
    variable = 'potLi'
    boundary = SE_left
    sort_by = y
  []
  [SESide_potential]
    type = SideValueSampler
    variable = 'potLi'
    boundary = SE_bottom
    sort_by = x
  []
  [SEInt_potential]
    type = SideValueSampler
    variable = 'potLi'
    boundary = interface0
    sort_by = x
  []
#  [SEInt_potential]
#    type = LineValueSampler
#    variable = 'potLi'
#    start_point = '0 0.5 0'
#    end_point = '0.5 0.5 0'
#    num_points = 200
#    sort_by = x
#  []
  [SEAnode_current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y'
    boundary = SE_left
    sort_by = y
  []
  [SEInt_current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y'
    boundary = interface0
    sort_by = x
  []

  [LiCathode_potential]
    type = SideValueSampler
    variable = 'potEn'
    boundary = Li_right
    sort_by = y
  []
  [LiSide_potential]
    type = SideValueSampler
    variable = 'potEn'
    boundary = Li_top
    sort_by = x
  []
#  [LiInt_potential]
#    type = SideValueSampler
#    variable = 'potEn'
#    boundary = interface0
#    block = blockLi
#    sort_by = x
#  []
  [LiInt_potential]
    type = LineValueSampler
    variable = 'potEn'
    start_point = '0  0.500001 0'
    end_point = '0.02 0.500001 0'
    num_points = 200
    sort_by = x
  []
  [LiAnode_current]
    type = SideValueSampler
    variable = 'iEn_x iEn_y'
    boundary = Li_left
    sort_by = y
  []

  [interfaceSE_normals]
    type = SideValueSampler
    variable = 'nodal_normal_x nodal_normal_y'
    boundary = interface0
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
    ionic_conductivity = 0.1
    metal_conductivity = 100000
    LiPotAnode = 0
    LiPotCathode = 0
  []
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    solve_type = 'NEWTON'
#  solve_type = 'PJFNK'
  [../]
[]

[Executioner]
  type = Steady

  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-8
  l_tol = 1e-4
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
  file_base = rst/mdlL2
  csv = true
[]
