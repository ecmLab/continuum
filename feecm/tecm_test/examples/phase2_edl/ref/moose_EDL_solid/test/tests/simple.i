e = 1.602e-19
kT = 4.142e-21
kTbye = 0.025855

[Mesh]
  active = 'mesh'
  [./mesh]
   dim = 1
   type = GeneratedMeshGenerator
   nx = 1e4
   xmax = 100e-6 # 100 micron thickness of solid electrolyte
  [../]
[]

[Variables]
  [./phi]
  order = FIRST
  family = LAGRANGE
  # initial_condition = 0.01
  [../]
  [./cpos]
  order	= FIRST
  family = LAGRANGE
 initial_condition = 1.1e-3
  [../]
  [./cneg]
  order	= FIRST
  family = LAGRANGE
  initial_condition = 1e-3
  [../]
[]

# [ICs]
  # [./phi]
  #   type = FunctionIC
  #   variable = phi
  #   min = 0.1
  #   max = 0.2
  # [../]
#   [./cpos]
#     type = RandomIC
#     variable = cpos
#     min = 1.0e-10
#     max = 1.1e-6
#   [../]
#   [./cneg]
#     type = RandomIC
#     variable = cneg
#     min = 1.0e-10
#     max = 1.1e-6
#   [../]
# []
  
[AuxVariables]
  [bounds_dummy]
    order = FIRST
    family = LAGRANGE
  []
[]
  
[Bounds]
  [cpos_upper_bound]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = cpos
    bound_type = upper
    bound_value = 1
  []
  [cpos_lower_bound]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = cpos
    bound_type = lower
    bound_value = 0
  []
  [cneg_upper_bound]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = cneg
    bound_type = upper
    bound_value = 1
  []
  [cneg_lower_bound]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = cneg
    bound_type = lower
    bound_value = 0
  []
[]

[Kernels]
  active = 'diff rhobyeps_term chem_pot_eqbm chempotvalue'
  [./diff]
    type = Diffusion
    variable = phi
  [../]
  [./rhobyeps_term]
    type = MaskedBodyForce
    variable = phi
    value = 1.0
    mask = rhobyeps
    coupled_variables = 'cpos cneg'
  [../]
    #will set mupos + muneg = 0
  [./chem_pot_eqbm]
    type = ChemPotEqual
    variable = cneg
    mua = mu_neg
    mub = mu_pos
  [../]
  #will set mu equal to mu_target
  [./chempotvalue]
    type = ChemPotValue
    variable = cpos
    mutarget = mutarget
    mu = mu_pos
  [../]
[]

[Materials]
  [./consts]
    type = ADGenericConstantMaterial
    prop_names =      'eps       N         zpos zneg  mu0_pos   mu0_neg    mutarget'
    prop_values = '1.328e-10 5.98412e28     1    1    0.06     0.065       0.0'
    #units         SI units    /m3         no    no     eV       eV          eV
  [../]
  [./chempot_pos] # in eV
     type = ADDerivativeParsedMaterial
     property_name = 'mu_pos'
     coupled_variables = 'phi cpos'
     material_property_names = 'mu0_pos zpos'
     expression = 'mu0_pos + ${kT}/${e} * log(cpos / (1-cpos)) + zpos * phi'
     derivative_order = 2
    outputs = exodus
  [../]
  [./chempot_neg] # in eV
     type = ADDerivativeParsedMaterial
     property_name = 'mu_neg'
     coupled_variables = 'phi cneg'
     material_property_names = 'mu0_neg zneg'
     expression = 'mu0_neg + ${kTbye} *  log(cneg / (1-cneg)) - zneg * phi'
     derivative_order = 2
     outputs = exodus
  [../]
  #Need to convert e and eps to non-AD for the next material
  [./convert_to_AD]
     type = MaterialADConverter
     ad_props_in = 'eps N'
     reg_props_out = 'eps_reg N_reg'
  [../]
  [./rhobyeps_f]
    type = DerivativeParsedMaterial
    property_name = 'rhobyeps'
    coupled_variables = 'phi cpos cneg'
    material_property_names = 'eps_reg N_reg'
    expression = '${e} * N_reg * ( cpos - cneg )/eps_reg'
    derivative_order = 2
    outputs = exodus
  [../]
[]

[BCs]
active = 'left rightn'
  [./left]
    type = DirichletBC
    variable = phi
    boundary = left
    value = -0.1
  [../]
  [./right]
    type = DirichletBC
    variable = phi
    boundary = 'right'
    value = 0.1
  [../]
  [./rightn]
    type = NeumannBC
    variable = phi
    boundary = right
    value = 0.0 #right side will have CNFL - bulk electroneutrality
  [../]
[]

[Adaptivity]
  marker = errorfrac # this specifies which marker from 'Markers' subsection to use
  steps = 20 # run adaptivity 2 times, recomputing solution, indicators, and markers each time
  # Use an indicator to compute an error-estimate for each element:
  [./Indicators]
    # create an indicator computing an error metric for the convected variable
    [./error]
      # arbitrary, use-chosen name
      type = GradientJumpIndicator
      variable = phi
      outputs = none
    [../]
  [../]
  # Create a marker that determines which elements to refine/coarsen based on error estimates
  # from an indicator:
  [./Markers]
    [./errorfrac]
      # arbitrary, use-chosen name (must match 'marker=...' name above
      type = ErrorFractionMarker
      indicator = error # use the 'error' indicator specified above
      refine = 0.2 # split/refine elements in the upper half of the indicator error range
      coarsen = 0 # don't do any coarsening
      outputs = none
    [../]
  [../]
[]

[Preconditioning]
  [./SMP]
  type = SMP
  full = true
#  petsc_options = '-snes_monitor -ksp_monitor_true_residual -snes_converged_reason -ksp_converged_reason'
  petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm      121                  preonly       lu           4'
  [../]
[]

[Executioner]
  automatic_scaling = true
  type = Steady
#  end_time = 10.0
#  dt = 0.1
#  scheme = bdf2
  verbose = True
  solve_type = 'Newton'
  l_max_its = 50
  l_tol = 1e-6
  nl_max_its = 50
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-14
  petsc_options = '-snes_monitor -ksp_monitor_true_residual -snes_converged_reason -ksp_converged_reason'
  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap -snes_type'
  petsc_options_value = 'asm      31                  preonly      lu          4 vinewtonssls'
#  line_search = 'none'
[]
			       
[Outputs]
  execute_on = 'final'
  exodus = true
[]

[Debug]
  show_material_props = true
  show_var_residual_norms = true
[]

