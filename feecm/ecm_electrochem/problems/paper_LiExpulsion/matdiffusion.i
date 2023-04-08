[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
  xmax = 1.0
  ymax = 1.0
  elem_type = QUAD4
[]

[Variables]
  [./u]
  [../]
[]

[Kernels]
  [./cres]
#    type = MatDiffusion
    type = Diffusion
    variable = u
#    diffusivity = Du
  [../]
  [./ctime]
    type = TimeDerivative
    variable = u
  [../]
[]

[Materials]
  [./Dc]
    type = DerivativeParsedMaterial
    f_name = Du
#    function = '0.01+u^2'
    function = '0.01'
    args = 'u'
    derivative_order = 1
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  [../]

  [./right]
    type = NeumannBC
    variable = u
    boundary = right
    value = -10
  [../]
[]

[VectorPostprocessors]
  [top_u]
    type = SideValueSampler
    variable = 'u'
    boundary = top
    sort_by = x
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  scheme = 'BDF2'
  dt = 1
  num_steps = 10
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
  csv = true
  file_base = rst/
[]
