
current_density = 10e-3
k_anode = 1e-3
faraday = 96.4853329
temperature = 298
gas_constant = 8.31446
c_ref = 0
c_ref_pore = 2.5925e-2
k_pore = 1e-3
c_init_Li = 7.69e-2
c_init_C = 1.6e-3
c_max_C = 3.05e-2
pressure = 1e-4
end_time = 150
[Mesh]
  [./mesh]
    type = FileMeshGenerator
    file = 'data/2blocks_50.e'
  [../]
  [./primary_subdomain]
    type = LowerDBlockFromSidesetGenerator
    input = mesh
    sidesets = 'right_left'
    new_block_name = 'left_right_primary_subdomain'
  [../]
  [./secondary_subdomain]
    type = LowerDBlockFromSidesetGenerator
    input = primary_subdomain
    sidesets = 'left_right'
    new_block_name = 'left_right_secondary_subdomain'
  [../]
[]
[Variables]
  [./li_conc]
    block = 'left1 right1'
  [../]
  [./equil]
    block = 'left_right_secondary_subdomain'
  [../]

[]

# [ICs]
#   [./li_conc_left]
#     type = FunctionIC
#     function = 'r:=sqrt((x-0.4)^2+(y-0.5)^2);if(r<0.05,5,1)'
#     variable = li_conc
#     block = 'left1'
#   [../]
#   [./li_conc_right]
#     type = ConstantIC
#     value = 0.8
#     variable = li_conc
#     block = 'right1'
#   [../]
# []


[Constraints]
  [./left_right_conc]
    type = EqualNormalFluxConstraint
    primary_boundary = 'right_left'
    secondary_boundary = 'left_right'
    primary_subdomain = 'left_right_primary_subdomain'
    secondary_subdomain = 'left_right_secondary_subdomain'
    primary_variable = 'li_conc'
    secondary_variable = 'li_conc'
    primary_tensor = false
    secondary_tensor = false
    variable = equil
    extra_vector_tags = 'ref'
    primary_mat_prop = diffusivity
    secondary_mat_prop = diffusivity
    use_displaced_mesh = false
    displacement = 'ux uy'
  [../]
[]

[Kernels]
  [./li_conc]
    type = ADChemoMechanoDiffusion
    variable = li_conc
    diffusivity = diffusivity
    use_displaced_mesh = false
    block ='left1 right1'
  [../]

  [./li_metal_dt]
    type = ADTimeDerivative
    variable = li_conc
    use_displaced_mesh = false
    block = 'left1 right1'
  [../]
[]

[Materials]
  [./diffusivity]
    type = ADIsotropicDiffusionMaterial
    diffusion_coef = 1
    block = 'left1'
  [../]
  [./diff]
    type = ADIsotropicDiffusionMaterial
    diffusion_coef = 100
    block = 'right1'
  [../]
[]

[BCs]
  [./conc_li_flux]
    type = ADNeumannBC
    boundary = 'left_left'
    variable = li_conc
    value = '1'
    extra_vector_tags = 'ref'
  [../]
  [./diffusion_flux_bc]
    type = MaterialDiffusionFluxBC
    boundary = 'left_right'
    variable = li_conc
    extra_vector_tags = 'ref'
  [../]
  [./diffusion_flux_bc1]
    type = MaterialDiffusionFluxBC
    boundary = 'right_left'
    variable = li_conc
    extra_vector_tags = 'ref'
  [../]

[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  automatic_scaling = true
  # compute_scaling_once = false
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1     NONZERO               1e-20               '
  dt = 1
  # num_steps = 20
  l_max_its = 50
  nl_max_its = 25
  # nl_abs_tol = 1e-4
  nl_rel_tol = 1e-6
  # end_time =
  num_steps = 10
  # scaling_group_variables = 'V flux_interlayer flux_pore'
  resid_vs_jac_scaling_param = 0.5
[]

[Problem]
  type = ReferenceResidualProblem
  # solution_variables = 'ux uy normal_lm thermal_lm V thermal_lm2 li_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'li_conc equil'
  acceptable_iterations = 2
  # coord_type = RZ
[]

[Postprocessors]
  [./equilibrating_flux]
    type = ElementAverageValue
    block = 'left_right_secondary_subdomain'
    variable = equil
  [../]
[]

[Outputs]
  exodus = true
  # csv = true
  [./csv]
    type = CSV
    file_base = csv/test_equil
  [../]
  [./out]
    type = Exodus
    file_base = rst/test_equil
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
[]
[AuxVariables]
  [./ux]
    block = 'left1 right1'
  [../]
  [./uy]
    block = 'left1 right1'
  [../]

  [./li_metal_flux_x]
    order = CONSTANT
    family = MONOMIAL
    block = 'right1'
  [../]
  [./li_metal_flux_y]
    order = CONSTANT
    family = MONOMIAL
    block = 'right1'
  [../]
  [./li_metal_flux_z]
    order = CONSTANT
    family = MONOMIAL
    block = 'right1'
  [../]
[]
[AuxKernels]
  [./ux]
    type = ConstantAux
    value = 0
    variable = ux
    block = 'left1 right1'
  [../]
  [./uy]
    type = ConstantAux
    value = 0
    variable = uy
    block = 'left1 right1'
  [../]

  [./li_metal_flux_x1]
    type = ADDiffusionFluxAux
    variable = li_metal_flux_x
    component = x
    diffusion_variable = li_conc
    block = 'right1'
    diffusivity = diffusivity
  [../]

  [./li_metal_flux_y1]
    type = ADDiffusionFluxAux
    variable = li_metal_flux_y
    component = y
    diffusion_variable = li_conc
    block = 'right1'
    diffusivity = diffusivity
  [../]
  [./li_metal_flux_z1]
    type = ADDiffusionFluxAux
    variable = li_metal_flux_z
    component = z
    diffusion_variable = li_conc
    block = 'right1'
    diffusivity = diffusivity
  [../]
[]
