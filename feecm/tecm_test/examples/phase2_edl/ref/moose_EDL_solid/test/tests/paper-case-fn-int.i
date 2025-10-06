echarge = 4.140E+01
kTbye = 1
V_scale = 2.58E-02
# kT = ${fparse e*kTbye}
dfepos = ${fparse 0.5/V_scale}  # 1.935E+01
Bpos = ${fparse 0./V_scale} #-7.739E+00 # -1.548E+01
etapos = 3.000E-03 # 0.0003
dfeneg = ${fparse dfepos} # 37.53478261 # in J/e = eV
Bneg = ${fparse 0./V_scale}
etaneg = 3.000E-03 #  0.0003
N = 5.000E+07
epsr = 10.0
eps0 = 5.913E+00
V0 =  ${fparse 1.0/V_scale} #1.935E+01 
L = 3.800E+01
c0 = ${fparse exp(-(dfepos+dfeneg)/2/kTbye)}
alpha = 0.1
fp = ${fparse 2.0/V_scale} #65.
fm = ${fparse 2.0/V_scale}  #65.
fpm = 0.0
Ld = ${fparse sqrt(epsr*eps0*kTbye/2/echarge/N/c0)}

[Mesh]
active = 'fmesh boundaries'
  [./mesh]
   dim = 1
   type = GeneratedMeshGenerator
   nx = 1e5
   xmax = ${L} # thickness of solid electrolyte
  [../]
  [./fmesh]
    type = FileMeshGenerator
    file = '/Users/zeeshan/projects/edl_solid/test/tests/line-1e5-L38.msh'
  [../]
  [./boundaries]
    type = SideSetsFromNormalsGenerator
    input = fmesh
    normals = '-1 0 0
                1 0 0'
    fixed_normal = true
    new_boundary = 'left right'
  [../]
[]

#For V=0.5, 1.0: 1.0e-2, 1.0, 1.0e4
[Variables]
  [./phi]
  order = FIRST
  family = LAGRANGE
  # initial_condition = 0 #${fparse -V0/2}
  # scaling = ${fparse 1.0e4 * exp(-V0)}
  # scaling = 1.0e-2
  [../]
  [./cpos]
  order	= FIRST
  family = LAGRANGE
  # initial_condition = 1.0 # ${fparse exp(-dfepos/kTbye)}
  # scaling = ${fparse exp(dfepos/kTbye)}
  # scaling = ${fparse 1e5 * exp(-V0)}
  # scaling = 1.0
  [../]
  [./cneg]
  order	= FIRST
  family = LAGRANGE
  # initial_condition = 1.0 # ${fparse exp(-dfeneg/kTbye)}
  # scaling = ${fparse exp(dfepos/kTbye)}
  # scaling = ${fparse 1e5 * exp(V0)}
  # scaling = 1.0e4
  [../]
[]

[ICs]
  [./phi]
    type = FunctionIC
    variable = phi
    function = phi_exp
  [../]
 [./cpos]
   type = FunctionIC
    variable = cpos
    function = cpos_linear
   [../]
   [./cneg]
     type = FunctionIC
     variable = cneg
     function = cneg_linear
   [../]
[]

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
    bound_value = ${fparse 1/c0}
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
    bound_value = ${fparse 1/c0}
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
  active = 'diff rhobyeps_term chempotvalue chempotvalue2'
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
  [./phi_linear]
    type = ParsedFunction
    expression = '${fparse -V0} * (1 - x/${L})'
  [../]
  [./phi_exp]
    type = ParsedFunction
    # expression = '${fparse -V0} * exp(-x/${Ld})'
    expression = '4 * ${kTbye} * atanh(tanh(-${V0}/4/${kTbye}) * exp(-x/${Ld}))'
  [../]
  [./cpos_linear]
    type = ParsedFunction
    # expression = '1/${c0} / (1 + 1/${c0} * exp(-${V0} * exp(-x/${Ld}) / ${kTbye}))'
    expression = '1/${c0} / (1/0.5 + 1/${c0} * exp( 4 * ${kTbye} * atanh(tanh(-${V0}/4/${kTbye}) * exp(-x/${Ld})) / ${kTbye}))'
    # expression = 'min(1.0 * exp(${V0} * exp(-x/${Ld}) / ${kTbye}), 0.8/${c0})'
    # expression = '1 + (min(${fparse exp(V0/kTbye)},${fparse 1/c0}) - 1) * exp(-x/${Ld})'
    #'${fparse exp(V0/kTbye)} + (1 - ${fparse exp(V0/kTbye)}) * x/${L}'
  [../]
  [./cneg_linear]
    type = ParsedFunction
    expression = '1/${c0} / (1/0.5 + 1/${c0} * exp( 4 * ${kTbye} * atanh(tanh(${V0}/4/${kTbye}) * exp(-x/${Ld})) / ${kTbye} ))'
#    expression = '1.0 * exp(-${V0} * exp(-x/${Ld}) / ${kTbye})'
    # expression = '1 + (${fparse exp(-V0/kTbye)} - 1) * exp(-x/${Ld})'
    #'${fparse exp(-V0/kTbye)} + (1 - ${fparse exp(-V0/kTbye)}) * x/${L}'
  [../]
  [./Vfunc]
    type = ParsedFunction
    expression = '${fparse -V0}*t'
  [../]
  [./dfepos0_func]
    type = ParsedFunction
    symbol_names = 'Apos  Bpos  etapos'
    symbol_values = '${dfepos} ${Bpos} ${etapos}'  #eV and cpos in per nm
    expression = 'if(etapos=0, Apos, Apos + Bpos * exp(-x/etapos))'
  [../]
  [./dfeneg0_func]
    type = ParsedFunction
    symbol_names = 'Aneg   Bneg  etaneg'
    symbol_values = '${dfeneg}  ${Bneg}  ${etaneg}'
    expression = 'if(etaneg=0, Aneg, Aneg + Bneg * exp(-x/etaneg))'
  [../]
[]

[Materials]
  [./consts]
    type = ADGenericConstantMaterial
    prop_names =      'eps                    N       zpos  zneg    mutarget  fp   fm  fpm poisson_multiplier alpha Ld dfe c0'
    prop_values =   '${fparse epsr*eps0}    ${N}     1     1        0.0       ${fp}  ${fm}  ${fpm}   ${fparse N*echarge*c0/epsr/eps0} ${alpha} ${Ld} ${dfepos} ${c0}' # ${fparse e*N/epsr/eps0*c0}
    #units             SI units             /m3        no    no        eV    eV   eV   eV
    outputs = exodus
  [../]
  [./mu0_pos_f]
    type = ADGenericFunctionMaterial
    prop_names = mu0_pos
    prop_values = dfepos0_func
    outputs = exodus
  [../]
  [./mu0_neg_f]
    type = ADGenericFunctionMaterial
    prop_names = mu0_neg
    prop_values = dfeneg0_func
    outputs = exodus
  [../]
  [./fint_pos]
    type = ADDerivativeParsedMaterial
    property_name = 'muint_pos'
    coupled_variables = 'phi cpos cneg'
    material_property_names = 'fp fpm'
    # expression = 'fp * cpos * cpos * ${c0}^2 + fpm * cneg * ${c0}'
    # expression = 'fp * 1/(1-cpos*${c0}) + fpm * cneg'
    # expression = 'fp * 4 * (cpos * ${c0} * cpos * ${c0}) / (1 + exp(-100*(cpos * ${c0} - 0.5))) + fpm * cneg * ${c0}'
    expression = 'fp * (cpos * ${c0})^2 / (1-cpos*${c0})^2 + fpm * cneg * ${c0}'
    derivative_order = 2
    outputs = exodus
  [../]
  [./fint_neg]
    type = ADDerivativeParsedMaterial
    property_name = 'muint_neg'
    coupled_variables = 'phi cpos cneg'
    material_property_names = 'fm fpm'
    # expression = 'fm * cneg * cneg *${c0}^2 + fpm * cpos * ${c0}'
    # expression = 'fm * 1/(1-cneg*${c0}) + fpm * cpos'
    # expression = 'fm * 4 * (cneg * ${c0} * cneg * ${c0}) / (1 + exp(-100*(cneg * ${c0} - 0.5))) + fpm * cpos * ${c0}'
    expression = 'fm * (cneg * ${c0})^2 / (1-cneg*${c0})^2 + fpm * cpos * ${c0}'
    derivative_order = 2
    outputs = exodus
  [../] 
  [./chempot_pos] # in eV
     type = ADDerivativeParsedMaterial
     property_name = 'mu_pos'
     coupled_variables = 'phi cpos cneg'
     material_property_names = 'mu0_pos zpos muint_pos'
    #  expression = 'exp(mu0_pos + zpos * phi + fp * cpos * cpos + fpm * cneg) * (cpos * ${c0} / (1 - cpos * ${c0}))^${kTbye} - 1'
    expression = 'if(cpos*${c0} < 0, 1e10*(cpos - 1e8), mu0_pos + ${kTbye} * log(cpos * ${c0} / (1-cpos*${c0})) + zpos * phi + muint_pos)'
    #  expression = 'if(mu0_pos + zpos * phi + fp * cpos * cpos + fpm * cneg < -35., cpos * ${c0} - 1/(exp(mu0_pos + zpos * phi) + 1), if(mu0_pos + zpos * phi > 35., cpos, mu0_pos - ${kTbye} *  log(cpos / (1-cpos)) + zpos * phi + fp * cpos * cpos + fpm * cneg))'
    #  expression = 'mu0_pos + ${kTbye} * ( log(cpos) + ${fparse -(dfepos+dfeneg)/2/kTbye} ) + zpos * phi + fp * cpos * cpos + fpm * cneg'
    #  expression = 'mu0_pos - ${kTbye} *  log(1 / cpos / ${c0} - 1 / 0.1)  + zpos * phi'
    # expression = 'if(mu0_pos + zpos * phi < -35., cpos * ${c0} - 0.1, if(mu0_pos + zpos * phi > 35., cpos, mu0_pos - ${kTbye} *  log(1 / cpos / ${c0} - 1 / 0.1) + zpos * phi + fp * cpos * cpos + fpm * cneg))'
    # expression = 'if(mu0_pos + zpos * phi < -35., cpos * ${c0} , if(mu0_pos + zpos * phi > 35., cpos, mu0_pos - ${kTbye} *  log(cpos / (1-cpos)) + zpos * phi + fp * cpos * cpos + fpm * cneg))'
    #  derivative_order = 2
    outputs = exodus
  [../]
  [./chempot_neg] # in eV
     type = ADDerivativeParsedMaterial
     property_name = 'mu_neg'
     coupled_variables = 'phi cneg cpos'
     material_property_names = 'mu0_neg zneg muint_neg'
    #  expression = 'exp(mu0_neg - zneg * phi + fm * cneg * cneg + fpm * cpos) * (cneg * ${c0} / (1 - cneg * ${c0}))^${kTbye} - 1'
    expression = 'if(cneg*${c0} < 0, 1e10*(cneg - 1e-8), mu0_neg + ${kTbye} * log(cneg * ${c0} / (1-cneg*${c0})) - zneg * phi + muint_neg)'
    #  expression = 'mu0_neg + ${kTbye} * ( log(cneg) + ${fparse -(dfepos+dfeneg)/2/kTbye} )  - zneg * phi + fm * cneg * cneg + fpm * cpos'
    #  expression = 'mu0_neg + ${kTbye} * log(1/cneg/${c0} - 1/0.1) - zneg * phi'
    # expression = 'if(mu0_neg - zneg * phi < -35.0, cneg * ${c0} - 0.1, if(mu0_neg - zneg * phi > 35., cneg, mu0_neg - ${kTbye} * log(1 / cneg / ${c0} - 1 / 0.1) - zneg * phi + fm * cneg * cneg + fpm * cpos))'
    # expression = 'if(mu0_neg - zneg * phi < -35.0, cneg * ${c0} , if(mu0_neg - zneg * phi > 35., cneg, mu0_neg + ${kTbye} * log(cneg / (1-cneg)) - zneg * phi + fm * cneg * cneg + fpm * cpos))'
     derivative_order = 2
     outputs = exodus
  [../]
  #Need to convert e and eps to non-AD for the next material
  [./convert_to_AD]
     type = MaterialADConverter
     ad_props_in = 'eps N poisson_multiplier'
     reg_props_out = 'eps_reg N_reg poisson_multiplier_reg'
  [../]
  [./rhobyeps_f]
    type = DerivativeParsedMaterial
    property_name = 'rhobyeps'
    coupled_variables = 'phi cpos cneg'
    material_property_names = 'eps_reg N_reg poisson_multiplier_reg'
    # expression = '${e} * N_reg * ( cpos - cneg )/eps_reg * ${c0}'
    expression = 'poisson_multiplier_reg * (cpos - cneg)'
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
    value = ${fparse -V0}
    [../]
  [./leftf]
    type = FunctionDirichletBC
    variable = phi
    boundary = left
    function = Vfunc
  [../]
  [./right]
    type = DirichletBC
    variable = phi
    boundary = 'right'
    value = 0 #${fparse -V0}
  [../]
  [./rightn]
    type = NeumannBC
    variable = phi
    boundary = right
    value = 0 #right side will have CNFL
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

# [Adaptivity]
#   marker = errorfrac # this specifies which marker from 'Markers' subsection to use
#   steps = 4 # run adaptivity 2 times, recomputing solution, indicators, and markers each time
#   max_h_level = 4 # maximum number of times to refine
#   cycles_per_step = 1 # number of times to refine/coarsen each step
#   # Use an indicator to compute an error-estimate for each element:
#   recompute_markers_during_cycles = true
#   [./Indicators]
#     # create an indicator computing an error metric for the convected variable
#     [./error]
#       # arbitrary, use-chosen name
#       type = GradientJumpIndicator
#       variable = phi
#       outputs = none
#     [../]
#   [../]
#   # Create a marker that determines which elements to refine/coarsen based on error estimates
#   # from an indicator:
#   [./Markers]
#     [./errorfrac]
#       # arbitrary, use-chosen name (must match 'marker=...' name above
#       type = ErrorFractionMarker
#       indicator = error # use the 'error' indicator specified above
#       refine = 0.5 # split/refine elements in the upper half of the indicator error range
#       coarsen = 0 # don't do any coarsening
#       outputs = none
#     [../]
#   [../]
# []

[Preconditioning]
  [./SMP]
  type = SMP
  full = true
  petsc_options = '-snes_monitor -ksp_monitor_true_residual -snes_converged_reason -ksp_converged_reason'
  petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap -sub_pc_factor_shift_type -snes_type'
  petsc_options_value = 'asm      121                  preonly       lu           4 NONZERO vinewtonssls'
  # petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter -snes_type'
  # petsc_options_value = '201                hypre    boomeramg      4 vinewtonssls'
  # petsc_options_iname = '-pc_type -snes_type'
  # petsc_options_value = 'lu vinewtonssls'
  [../]
[]

[Executioner]
  automatic_scaling = true
  resid_vs_jac_scaling_param = 0.5
  type = Steady
#   end_time = 1.0
#   dt = 0.05
#   dtmin = 1.0e-3
#   scheme = bdf2
  verbose = True
  solve_type = 'Newton'
  l_max_its = 100
  l_tol = 1e-6
  nl_max_its = 150
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-13
  petsc_options = '-snes_monitor -ksp_monitor_true_residual -snes_converged_reason -ksp_converged_reason'
  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap -sub_pc_factor_shift_type -snes_type'
  petsc_options_value =  'asm      31                  preonly         lu             4       NONZERO vinewtonssls'
  # petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter -snes_type'
  # petsc_options_value = '201                hypre    boomeramg      4 vinewtonssls'
  # petsc_options_iname = '-pc_type -snes_type'
  # petsc_options_value = 'lu vinewtonssls'
#  line_search = 'none'
#   [TimeStepper]
#     type = IterationAdaptiveDT
#     dt = 0.05
#     growth_factor = 2
#   []
[]

[Outputs]
  execute_on = 'TIMESTEP_END' #'final' 'initial failed' 
  exodus = true
[]

[Debug]
  # show_material_props = true
  show_var_residual_norms = true
[]
