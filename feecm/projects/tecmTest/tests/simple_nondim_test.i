# Simple test of NonDimensionalParameters material
# Copyright 2025, CEWLAB, All Rights Reserved

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 2
[]

[Variables]
  [dummy]
  []
[]

[Materials]
  [nondim_params]
    type = NonDimensionalParameters
    c0 = 41528.0
    T0 = 298.0
    D0 = 2.5e-12
    j0 = 1.0
    Omega0 = 1.304e-5
    F = 96485.0
    R = 8.314
  []
[]

[Kernels]
  [dummy_diffusion]
    type = Diffusion
    variable = dummy
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = dummy
    boundary = left
    value = 1
  []
  [right]
    type = DirichletBC
    variable = dummy
    boundary = right
    value = 0
  []
[]

[Postprocessors]
  [D0_check]
    type = ElementAverageMaterialProperty
    mat_prop = D0
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
[]