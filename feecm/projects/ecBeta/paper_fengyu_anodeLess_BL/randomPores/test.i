[MeshGenerators]
 [GeneratedMesh]
   type = FileMeshGenerator
   file = data/mdlR3.msh
 []

 [./inter1]
   type = SideSetsBetweenSubdomainsGenerator
   input = GeneratedMesh
   master_block = 'blockSE'
   paired_block = 'blockLi1 blockLi2'
   new_boundary = 'interface01'
 [../]

# [./inter2]
#   type = SideSetsBetweenSubdomainsGenerator
#   input = GeneratedMesh
#   master_block = 'blockSE'
#   paired_block = 'blockLi2'
#   new_boundary = 'interface02'
# [../]

 [./break_boundary1]
   input = inter1
   type = BreakBoundaryOnSubdomainGenerator
 [../]

# [./break_boundary2]
#   input = inter2
#   type = BreakBoundaryOnSubdomainGenerator
# [../]

[]

[Mesh]
  type = MeshGeneratorMesh
[]

[Variables]
  [potLi]
   block = 'blockSE'
  []
  [potEn]
   block = 'blockSE'
  []
  [potMt]
   block = 'blockLi1 blockLi2'
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
    block = 'blockLi1 blockLi2'
  []
  [iMt_y]
    order = CONSTANT
    family = MONOMIAL
    block = 'blockLi1 blockLi2'
  []
  [potDif]
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
    block = 'blockLi1 blockLi2'
  []
[]

[InterfaceKernels]
  [interface1_current]
    type = InterfaceDiffusion
    variable = potEn
    neighbor_var = potMt
    potLi    = potLi
    boundary = interface01
  []
  [interface1_potential]
    type = PenaltyInterfaceDiffusion
    variable = potEn
    neighbor_var = potMt
    boundary = interface01
    penalty = 1e10
  []
#  [interface2]
#    type = InterfaceDiffusion
#    variable = potEn
#    neighbor_var = potMt
#    potLi    = potLi
#    boundary = interface02
#  []
#  [interface2_potential]
#    type = PenaltyInterfaceDiffusion
#    variable = potEn
#    neighbor_var = potMt
#    boundary = interface02
#    penalty = 1e10
#  []
[]

[AuxKernels]
  [iLi_x]
    type = CurrentDensityLi
    variable = iLi_x
    component = x
    execute_on = timestep_end
    potential = potLi
    block = blockSE
  []
  [iLi_y]
    type = CurrentDensityLi
    variable = iLi_y
    component = y
    execute_on = timestep_end
    potential = potLi
    block = blockSE
  []
  [iEn_x]
    type = CurrentDensityEn
    variable = iEn_x
    component = x
    execute_on = timestep_end
    potential = potEn
    block = blockSE
  []
  [iEn_y]
    type = CurrentDensityEn
    variable = iEn_y
    component = y
    execute_on = timestep_end
    potential = potEn
    block = blockSE
  []
  [iMt_x]
    type = CurrentDensityMt
    variable = iMt_x
    component = x
    execute_on = timestep_end
    potential = potMt
    block = 'blockLi1 blockLi2'
  []
  [iMt_y]
    type = CurrentDensityMt
    variable = iMt_y
    component = y
    execute_on = timestep_end
    potential = potMt
    block = 'blockLi1 blockLi2'
  []
  [potDif]
    type = OverPotential
    variable = potDif
    execute_on = timestep_end
    potLi = potLi
    potEn = potEn
  []
[]

[BCs]
  [inletLi_BV]
    type = DepositionLi
    variable = potLi
    potEn    = potEn
    boundary = bottom_to_blockSE
  []
  [outletLi_BV]
    type = DepositionLi
    variable = potLi
    potEn = potEn
    boundary = top_to_blockSE
  []
  [Interface1Li_BV]
    type = DepositionLi
    variable = potLi
    potEn = potEn
    boundary = interface1_to_blockSE
  []
  [Interface2Li_BV]
    type = DepositionLi
    variable = potLi
    potEn = potEn
    boundary = interface2_to_blockSE
  []
  [Interface3Li_BV]
    type = DepositionLi
    variable = potLi
    potEn = potEn
    boundary = interface3_to_blockSE
  []

  [outletEn_current]
    type = NeumannEn
    variable = potEn
    potLi    = potLi
    boundary = bottom_to_blockSE
  []
  [interface3En_current]
    type = DepositionEn
    variable = potEn
    potLi = potLi
    boundary = interface3_to_blockSE
  []
  [inletEn_potential]
    type = DirichletBC
    variable = potEn
    boundary = top_to_blockSE
    value = 0
  []

[]

[Postprocessors]
  [inlet_En]
    type = SideTotEn
    variable = potEn
    boundary = top_to_blockSE
    diffusivity = electronic_conductivity
  []
  [interface1SE_En]
    type = SideTotEn
    variable = potEn
    boundary = interface1_to_blockSE
    diffusivity = electronic_conductivity
  []
  [interface2SE_En]
    type = SideTotEn
    variable = potEn
    boundary = interface2_to_blockSE
    diffusivity = electronic_conductivity
  []
  [interface3SE_En]
    type = SideTotEn
    variable = potEn
    boundary = interface3_to_blockSE
    diffusivity = electronic_conductivity
  []
  [interface1Li_En]
    type = SideTotEn
    variable = potMt
    boundary = interface1_to_blockLi1
    diffusivity = metal_conductivity
  []
  [interface2Li_En]
    type = SideTotEn
    variable = potMt
    boundary = interface2_to_blockLi2
    diffusivity = metal_conductivity
  []
  [inner1_En]
    type = SideTotEn
    variable = potMt
    boundary = inner1_to_blockLi1
    diffusivity = metal_conductivity
  []
  [inner2_En]
    type = SideTotEn
    variable = potMt
    boundary = inner2_to_blockLi2
    diffusivity = metal_conductivity
  []
  [outlet_En]
    type = SideTotEn
    variable = potEn
    boundary = bottom_to_blockSE
    diffusivity = electronic_conductivity
  []
[]

[VectorPostprocessors]
  [top_potential]
    type = SideValueSampler
    variable = 'potLi potEn potMt'
    boundary = top_to_blockSE
    sort_by = y
  []
  [interface1SE_potential]
    type = SideValueSampler
    variable = 'potLi potEn potMt'
    boundary = interface1_to_blockSE
    sort_by = y
  []
  [interface2SE_potential]
    type = SideValueSampler
    variable = 'potLi potEn potMt'
    boundary = interface2_to_blockSE
    sort_by = y
  []
  [interface3SE_potential]
    type = SideValueSampler
    variable = 'potLi potEn potMt'
    boundary = interface3_to_blockSE
    sort_by = y
  []
  [interface1Li_potential]
    type = SideValueSampler
    variable = 'potLi potEn potMt'
    boundary = interface1_to_blockLi1
    sort_by = y
  []
  [interface2Li_potential]
    type = SideValueSampler
    variable = 'potLi potEn potMt'
    boundary = interface2_to_blockLi2
    sort_by = y
  []
  [bottom_potential]
    type = SideValueSampler
    variable = 'potLi potEn potMt'
    boundary = bottom_to_blockSE
    sort_by = y
  []

  [top_current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y iEn_x iEn_y iMt_x iMt_y'
    boundary = top_to_blockSE
    sort_by = y
  []
  [interface1SE_current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y iEn_x iEn_y iMt_x iMt_y'
    boundary = interface1_to_blockSE
    sort_by = y
  []
  [interface2SE_current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y iEn_x iEn_y iMt_x iMt_y'
    boundary = interface2_to_blockSE
    sort_by = y
  []
  [interface3SE_current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y iEn_x iEn_y iMt_x iMt_y'
    boundary = interface3_to_blockSE
    sort_by = y
  []
  [interface1Li_current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y iEn_x iEn_y iMt_x iMt_y'
    boundary = interface1_to_blockLi1
    sort_by = y
  []
  [interface2Li_current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y iEn_x iEn_y iMt_x iMt_y'
    boundary = interface2_to_blockLi2
    sort_by = y
  []
  [inner1_current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y iEn_x iEn_y iMt_x iMt_y'
    boundary = inner1_to_blockLi1
    sort_by = y
  []
  [inner2_current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y iEn_x iEn_y iMt_x iMt_y'
    boundary = inner2_to_blockLi2
    sort_by = y
  []
  [bottom_current]
    type = SideValueSampler
    variable = 'iLi_x iLi_y iEn_x iEn_y iMt_x iMt_y'
    boundary = bottom_to_blockSE
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
    inlet_current = 1
    exchange_current = 30
    reaction_rate = 0.5
    ionic_conductivity = 0.5
#    electronic_conductivity = 0.00001
    metal_conductivity = 10
    electronic_conductivity = 0.00001
#    metal_conductivity = 1
  []
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    solve_type = 'NEWTON'
  [../]
[]

[Executioner]
  type = Steady
#  solve_type = NEWTON

  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-8
#  nl_max_its = 6
  l_tol = 1e-4
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
  file_base = rst/mdlR3
  csv = true
[]
