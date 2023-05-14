name = 't1'

[Mesh]
 [GeneratedMesh]
   type = FileMeshGenerator
   file = data/${name}.msh
 []
 
[]

[Problem]
  coord_type = RZ
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
    type = ParamDiffusion
    variable = potLi
    conductivity = ionic_conductivity
#    block = blockSE
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
[]

[BCs]
  [Li_anode_potential]
    type = DirichletBC
    variable = potLi
    boundary = And_btm
    value = 0
  []
  [SE_anode_BV]
    type = SingleSEElectrodeBV
    variable = potLi
    LiPotElectrode = LiPotAnode
    boundary = And_btm
  []
  [SE_cathode_current]
    type = SingleSEElectrodeNeumann
    variable = potLi
    boundary = SE_btm
  []

[]


#[Postprocessors]
#  [SEAnd_LiTotCurrent]
#    type = SideTotCrnt
#  variable = potLi
#    boundary = 'And_btm'
#    conductivity = ionic_conductivity
#  []
#  [SECtd_LiTotCurrent]
#    type = SideTotCrnt
#    variable = potLi
#    boundary = SE_btm
#    conductivity = ionic_conductivity
#  []
#[]

[VectorPostprocessors]
  [SEAnd_cntLi]
    type = SideValueSampler
    variable = 'iLi_x iLi_y'
    boundary = 'And_btm'
    sort_by = x
  []
[]

[Materials]
## Unit system used in this code:
## length: um,     potential: V,    current: nA
## current density: mA/cm^2,  conductivity: mS/cm
## Convert all parameters to these units !!!

  [packed]
    type = SingleSE
    applied_current = 0.2
    exchange_current = 1.3
    reaction_rate = 0.5
    ionic_conductivity = 0.1
    metal_conductivity = 10000
    gb_conductivity = 0.001
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
  nl_max_its = 50
  l_tol = 1e-6
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = false
  file_base = rst/${name}
  csv = true
[]
