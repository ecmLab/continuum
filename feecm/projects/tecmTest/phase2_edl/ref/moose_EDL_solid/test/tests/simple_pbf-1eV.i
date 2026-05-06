e = 1.602e-19
kT = 4.142e-21
kTbye = 0.025855
dfepos = 1.0
dfeneg = 0.9
N = 5.98412e28
epsr = 15.0
eps0 = 8.854e-12
V0 = -0.1
L = 100e-4
alpha = 0.1

[Mesh]
  active = 'mesh'
  [./mesh]
   dim = 1
   type = GeneratedMeshGenerator
   nx = 1e5
   xmax = ${L} # 100 micron thickness of solid electrolyte
  [../]
[]

[Variables]
  [./phi]
  order = FIRST
  family = LAGRANGE
  initial_condition = ${fparse V0/2}
  [../]
  [./cpos]
  order	= FIRST
  family = LAGRANGE
  initial_condition = ${fparse exp(-dfepos/kTbye)}
  [../]
  [./cneg]
  order	= FIRST
  family = LAGRANGE
  initial_condition = ${fparse exp(-dfeneg/kTbye)}
  [../]
[]

# [ICs]
#   [./phi]
#     type = RandomIC
#     variable = phi
#     min = 0.1
#     max = 0.2
#   [../]
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
  [./grad_phi_x]
    order = FIRST
    family = MONOMIAL
    outputs = exodus
  [../]
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
    #will set mupos = -muneg
  [./chem_pot_eqbm]
    type = ChemPotEqual
    variable = cneg
    mua = mu_neg
    mub = mu_pos
  [../]
  #will set mu equal to mu_target
  [./chempotvalue2]
    type = ChemPotValue
    variable = cneg
    mutarget = mutarget
    mu = mu_neg
  [../]
  [./chempotvalue]
    type = ChemPotValue
    variable = cpos
    mutarget = mutarget
    mu = mu_pos
  [../]
[]

[Functions]
  [./dfepos0_func]
    type = ParsedFunction
    symbol_names = 'Apos   Bpos   Cpos'
    symbol_values = '${dfepos} 0.  0.0'
#    symbol_values = '0.3e-19 -0.25e-19 5.30635e+08'  # Apos, Bpos change to eV, cpos in per nm
    expression = 'Apos + Bpos * exp(-x*1.0e9*Cpos)'
  [../]
  [./dfeneg0_func]
    type = ParsedFunction
    symbol_names = 'Aneg   Bneg  Cneg'
    symbol_values = '${dfeneg}  0.     0.'
#    symbol_values = '0.3e-19 -0.25e-19 3.01e+09'
    expression = 'Aneg + Bneg * exp(-x*1.0e9*Cneg)'
  [../]
[]

[Materials]
  [./consts]
    type = ADGenericConstantMaterial
    prop_names =      'eps                    N       zpos  zneg    mutarget  fp   fm  fpm'
    prop_values =   '${fparse epsr*eps0}      ${N}     1     1        0.0     0.   0.   0.'
    #units             SI units             /m3        no    no        eV    eV   eV   eV
  [../]
  [./mu0_pos_f]
    type = ADGenericFunctionMaterial
    prop_names = mu0_pos
    prop_values = dfepos0_func
  [../]
  [./mu0_neg_f]
    type = ADGenericFunctionMaterial
    prop_names = mu0_neg
    prop_values = dfeneg0_func
  [../]
  [./chempot_pos] # in eV
     type = ADDerivativeParsedMaterial
     property_name = 'mu_pos'
     coupled_variables = 'phi cpos cneg'
     material_property_names = 'mu0_pos zpos fp fpm'
     expression = 'mu0_pos + ${kTbye} * log(cpos / ${alpha} / (1 -  cpos / ${alpha})) + zpos * phi + fp * cpos  + fpm * cneg '
     derivative_order = 2
    outputs = exodus
  [../]
  [./chempot_neg] # in eV
     type = ADDerivativeParsedMaterial
     property_name = 'mu_neg'
     coupled_variables = 'phi cneg cpos'
     material_property_names = 'mu0_neg zneg fm fpm'
     expression = 'mu0_neg + ${kT}/${e} *  log(cneg / ${alpha} / (1-cneg / ${alpha} )) - zneg * phi + fm * cneg  + fpm * cpos'
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
    value = ${V0}
  [../]
  [./right]
    type = DirichletBC
    variable = phi
    boundary = 'right'
    value = 0.0
  [../]
  [./rightn]
    type = NeumannBC
    variable = phi
    boundary = right
    value = 0.0 #right side will have CNFL
  [../]
[]

[AuxKernels]
   [./dphidx]
    type = VariableGradientComponent
    variable = grad_phi_x
    component = x
    gradient_variable = phi
   [../]
[]

[Adaptivity]
  marker = errorfrac # this specifies which marker from 'Markers' subsection to use
  steps = 25 # run adaptivity 2 times, recomputing solution, indicators, and markers each time
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
#   petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
#   petsc_options_value = 'asm      121                  preonly       lu           8'
    petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
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
  nl_abs_tol = 1e-12
  # petsc_options = '-snes_monitor -ksp_monitor_true_residual -snes_converged_reason -ksp_converged_reason'
  # petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap -snes_type'
  # petsc_options_value =   'asm      31                  preonly         lu             8      vinewtonssls'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
#  line_search = 'none'
[]

[Outputs]
  execute_on = 'final'
  exodus = true
[]

[Debug]
  # show_material_props = true
  show_var_residual_norms = true
[]

