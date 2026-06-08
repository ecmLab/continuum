# Pressure-free monolithic battery — preliminary contact-stability model.
#
# Geometry: 2D plane-strain cross-section through the layered pouch.
#   Bottom -> top:  pkg_bot | cu | anode | separator (GPE) | cathode | al | pkg_top
#
# Eigenstrains: polymerization shrink in the GPE (eps_GPE) and in the
# package layers (eps_pkg). External BCs: lateral rollers (plane-strain
# symmetry), bottom roller, prescribed top traction sigma_stack (>= 0
# = compressive; 0 = pressure-free).
#
# Outputs: contact pressure on each interface, contact resistance Rc(pc),
# contact-stability metric Pi_c = <pc> / pmin.
#
# Sweep at the command line, e.g.:
#   ../../../eel-opt -i pressure_free_stack.i \
#     eps_GPE=0.05 eps_pkg=0.05 sigma_stack=0.0

# ------------------------------------------------------------------
# Geometry [mm]
# Layered active stack is enveloped by a polymer package on all four sides.
# The side panels of the package, with lateral rollers, contract in y and
# pull the top/bottom strips inward, compressing the inner stack.
# ------------------------------------------------------------------
W_inner   = 1.0
t_pkg     = 0.100
t_cu      = 0.025
t_anode   = 0.055
t_sep     = 0.055
t_cathode = 0.137
t_al      = 0.025

# ------------------------------------------------------------------
# Polymerization shrinkages (volumetric, in [0,1])
# ------------------------------------------------------------------
eps_GPE = 0.05
eps_pkg = 0.05

# ------------------------------------------------------------------
# External stack pressure [MPa] (0 = pressure-free)
# ------------------------------------------------------------------
sigma_stack = 0.0

# ------------------------------------------------------------------
# Lithiation breathing
#   alpha_Li_anode  > 0: graphite expands ~10% on full lithiation (soc -> 1)
#   alpha_Li_cathode< 0: cathode contracts on lithiation (or expands on
#                        delithiation/charging). LVPF is low-volume change (~5%).
#   cycle_amplitude  : amplitude of soc oscillation after cure ends (0 = static)
#   t_cycle          : period of one charge/discharge cycle [pseudo-time]
# Cure ramps from 0 to 1 over [0, 1] then plateaus; soc cycles only afterwards.
# ------------------------------------------------------------------
alpha_Li_anode   = 0.10
alpha_Li_cathode = -0.05
cycle_amplitude  = 0.0
t_cycle          = 1.0
end_time         = 1.0

# ------------------------------------------------------------------
# Operating temperature (applied post-cure)
#   Cure is treated as isothermal at T_ref (room temperature), so all
#   prior chapters remain unaffected.  After cure (t > 1), T ramps
#   linearly from T_ref to T_op over a half-unit transition window and
#   then holds at T_op for the remainder of the run.  T_op is the
#   ambient temperature at which the battery operates; sweeping T_op
#   over the DARPA Promethean envelope (-20..+60 C) at sigma_stack = 0
#   produces the temperature-dependent design map.  Setting T_op = T_ref
#   recovers the isothermal baseline.
#
# CTEs are volumetric (= 3 x linear CTE for isotropic homogenization).
# Linear-CTE references [ppm/K]:
#   Cu 17, Al 23, graphite-binder anode 5, GPE separator 150,
#   NMC-binder cathode 13, polymer + ceramic package 60.
# ------------------------------------------------------------------
T_ref           = 25.0     # [C] reference (room / cure) temperature
T_op            = 25.0     # [C] post-cure operating temperature
t_T_ramp        = 0.5      # [pseudo-time] cure-end to T_op transition window
alpha_T_cu      = 5.1e-5   # 3 * 17 ppm/K
alpha_T_anode   = 1.5e-5   # 3 * 5  ppm/K
alpha_T_sep     = 4.5e-4   # 3 * 150 ppm/K
alpha_T_cathode = 3.9e-5   # 3 * 13 ppm/K
alpha_T_al      = 6.9e-5   # 3 * 23 ppm/K
alpha_T_pkg     = 1.8e-4   # 3 * 60 ppm/K

# ------------------------------------------------------------------
# Moduli [MPa] and Poisson ratios
# ------------------------------------------------------------------
E_cu = 5.0e4
nu_cu = 0.34
E_anode = 1.0e3
nu_anode = 0.30
E_sep = 10.0
nu_sep = 0.45
E_cathode = 5.0e3
nu_cathode = 0.30
E_al = 3.0e4
nu_al = 0.33
E_pkg = 1.0e3
nu_pkg = 0.40

# ------------------------------------------------------------------
# Contact-resistance closure: Rc = R0 * exp(-alpha * pc)
# Pi_c = <pc> / pmin, with pmin chosen so Rc(pmin) = Rc_target.
#
# Defaults reflect literature for binder-bonded battery interfaces:
#   R0    = 200 mOhm.cm^2 (uncompressed PTFE-binder-bonded contact)
#   alpha = 0.5 /MPa      (pressure sensitivity, mid of 0.3-1.0 /MPa range)
#   Rc_target = 20 mOhm.cm^2 (proposal target, S1.5)
# These should be re-fit from the team's peel + 4-point measurements.
# ------------------------------------------------------------------
R0_mOhm_cm2 = 200.0
alpha_per_MPa = 0.5
Rc_target_mOhm_cm2 = 20.0
pmin = ${fparse log(R0_mOhm_cm2 / Rc_target_mOhm_cm2) / alpha_per_MPa}

[GlobalParams]
  energy_densities = 'dot(psi_m)'
  deformation_gradient = F
  mechanical_deformation_gradient = Fm
  eigen_deformation_gradient = Fg
[]

[Problem]
  kernel_coverage_check = false
  material_coverage_check = false
[]

[Mesh]
  [stack]
    type = CartesianMeshGenerator
    dim = 2
    dx = '${t_pkg} ${W_inner} ${t_pkg}'
    ix = '2 4 2'
    dy = '${t_pkg} ${t_cu} ${t_anode} ${t_sep} ${t_cathode} ${t_al} ${t_pkg}'
    iy = '2 2 6 2 6 2 2'
    # Row-major (bottom-to-top); 0 = package envelope, 1..5 = active stack layers
    subdomain_id = '
      0 0 0
      0 1 0
      0 2 0
      0 3 0
      0 4 0
      0 5 0
      0 0 0'
  []
  [rename]
    type = RenameBlockGenerator
    input = stack
    old_block = '0 1 2 3 4 5'
    new_block = 'package cu anode separator cathode al'
  []
  [if_cu_anode]
    type = SideSetsBetweenSubdomainsGenerator
    input = rename
    primary_block = 'cu'
    paired_block = 'anode'
    new_boundary = 'if_cu_anode'
  []
  [if_anode_sep]
    type = SideSetsBetweenSubdomainsGenerator
    input = if_cu_anode
    primary_block = 'anode'
    paired_block = 'separator'
    new_boundary = 'if_anode_sep'
  []
  [if_sep_cathode]
    type = SideSetsBetweenSubdomainsGenerator
    input = if_anode_sep
    primary_block = 'separator'
    paired_block = 'cathode'
    new_boundary = 'if_sep_cathode'
  []
  [if_cathode_al]
    type = SideSetsBetweenSubdomainsGenerator
    input = if_sep_cathode
    primary_block = 'cathode'
    paired_block = 'al'
    new_boundary = 'if_cathode_al'
  []
  use_displaced_mesh = false
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]

[AuxVariables]
  [cure]
    [AuxKernel]
      type = FunctionAux
      function = cure_func
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
  []
  [soc]
    [AuxKernel]
      type = FunctionAux
      function = soc_func
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
  []
  [T]
    [AuxKernel]
      type = FunctionAux
      function = T_func
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
  []
  [T_ref_var]
    initial_condition = ${T_ref}
  []
  [sigma_yy]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = ADRankTwoAux
      rank_two_tensor = cauchy
      index_i = 1
      index_j = 1
      execute_on = 'INITIAL TIMESTEP_END'
    []
  []
  [sigma_xx]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = ADRankTwoAux
      rank_two_tensor = cauchy
      index_i = 0
      index_j = 0
      execute_on = 'INITIAL TIMESTEP_END'
    []
  []
[]

[Functions]
  # Cure ramps from 0 to 1 over [0, 1], then plateaus.
  [cure_func]
    type = PiecewiseLinear
    x = '0 1 1.0e6'
    y = '0 1 1'
  []
  # State of charge: 0 during cure, then sinusoidal cycling between 0 and
  # cycle_amplitude with period t_cycle once cure is complete (t >= 1).
  [soc_func]
    type = ParsedFunction
    expression = 'if(t < 1.0, 0.0, ${cycle_amplitude} * 0.5 * (1 - cos(2*pi*(t-1.0)/${t_cycle})))'
  []
  # Temperature profile: T = T_ref during cure (t <= 1), then a linear
  # ramp from T_ref to T_op over t in [1, 1 + t_T_ramp], then held at
  # T_op for the remainder of the run.  T(1) = T_ref is enforced exactly,
  # so all cure-era physics is unchanged.
  [T_func]
    type = ParsedFunction
    expression = 'if(t < 1.0, ${T_ref}, ${T_ref} + (${T_op}-${T_ref})*min(1.0, (t-1.0)/${t_T_ramp}))'
  []
  # Top traction in y. Compressive (downward) is negative in this convention;
  # the function ramps with cure and then holds at the prescribed level.
  [sigma_stack_func]
    type = PiecewiseLinear
    x = '0 1 1.0e6'
    y = '0 ${fparse -1.0 * sigma_stack} ${fparse -1.0 * sigma_stack}'
  []
[]

[Kernels]
  [momentum_x]
    type = RankTwoDivergence
    variable = disp_x
    component = 0
    tensor = pk1
    factor = -1
  []
  [momentum_y]
    type = RankTwoDivergence
    variable = disp_y
    component = 1
    tensor = pk1
    factor = -1
  []
[]

[BCs]
  [fix_x_left]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  []
  [fix_x_right]
    type = DirichletBC
    variable = disp_x
    boundary = 'right'
    value = 0
  []
  [fix_y_bot]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  []
  [stack_pressure]
    type = ADFunctionNeumannBC
    variable = disp_y
    boundary = 'top'
    function = sigma_stack_func
  []
[]

[Materials]
  # Per-block elastic constants.
  [lambda]
    type = ADPiecewiseConstantByBlockMaterial
    prop_name = 'lambda'
    subdomain_to_prop_value = '
      package   ${fparse E_pkg*nu_pkg/(1+nu_pkg)/(1-2*nu_pkg)}
      cu        ${fparse E_cu*nu_cu/(1+nu_cu)/(1-2*nu_cu)}
      anode     ${fparse E_anode*nu_anode/(1+nu_anode)/(1-2*nu_anode)}
      separator ${fparse E_sep*nu_sep/(1+nu_sep)/(1-2*nu_sep)}
      cathode   ${fparse E_cathode*nu_cathode/(1+nu_cathode)/(1-2*nu_cathode)}
      al        ${fparse E_al*nu_al/(1+nu_al)/(1-2*nu_al)}'
  []
  [G]
    type = ADPiecewiseConstantByBlockMaterial
    prop_name = 'G'
    subdomain_to_prop_value = '
      package   ${fparse E_pkg/2/(1+nu_pkg)}
      cu        ${fparse E_cu/2/(1+nu_cu)}
      anode     ${fparse E_anode/2/(1+nu_anode)}
      separator ${fparse E_sep/2/(1+nu_sep)}
      cathode   ${fparse E_cathode/2/(1+nu_cathode)}
      al        ${fparse E_al/2/(1+nu_al)}'
  []
  # Per-block volumetric shrinkage (zero outside GPE and package).
  [eps_sh]
    type = ADPiecewiseConstantByBlockMaterial
    prop_name = 'eps_sh_vol'
    subdomain_to_prop_value = '
      package   ${eps_pkg}
      cu        0
      anode     0
      separator ${eps_GPE}
      cathode   0
      al        0'
  []
  # Per-block lithiation expansion coefficient (zero outside electrodes).
  [alpha_Li]
    type = ADPiecewiseConstantByBlockMaterial
    prop_name = 'alpha_Li'
    subdomain_to_prop_value = '
      package   0
      cu        0
      anode     ${alpha_Li_anode}
      separator 0
      cathode   ${alpha_Li_cathode}
      al        0'
  []
  # Per-block volumetric coefficient of thermal expansion [1/K].
  [alpha_T]
    type = ADPiecewiseConstantByBlockMaterial
    prop_name = 'alpha_T'
    subdomain_to_prop_value = '
      package   ${alpha_T_pkg}
      cu        ${alpha_T_cu}
      anode     ${alpha_T_anode}
      separator ${alpha_T_sep}
      cathode   ${alpha_T_cathode}
      al        ${alpha_T_al}'
  []

  # Kinematics: compose two eigen deformation gradients (cure shrink * Li breathing),
  # then pull back to the mechanical deformation gradient Fm.
  [Fc]
    type = CureShrinkageDeformationGradient
    cure_shrinkage_deformation_gradient = Fc
    cure = cure
    volumetric_shrinkage = eps_sh_vol
  []
  [Fl]
    type = LithiationDeformationGradient
    lithiation_deformation_gradient = Fl
    state_of_charge = soc
    lithiation_expansion_coefficient = alpha_Li
  []
  [Ft]
    type = ThermalDeformationGradient
    thermal_deformation_gradient = Ft
    temperature = T
    reference_temperature = T_ref_var
    CTE = alpha_T
  []
  [Fg_eigen]
    type = CompositeDeformationGradient
    composite_deformation_gradient = Fg_eigen
    factors = 'Fc Fl Ft'
  []
  [Fm]
    type = MechanicalDeformationGradient
    displacements = 'disp_x disp_y'
    swelling_deformation_gradient = Fg_eigen
  []

  # Constitutive: Neo-Hookean -> PK1 -> Cauchy.
  [neohookean]
    type = NeoHookeanSolid
    elastic_energy_density = psi_m
    lambda = lambda
    shear_modulus = G
  []
  [pk1]
    type = FirstPiolaKirchhoffStress
    first_piola_kirchhoff_stress = pk1
    deformation_gradient_rate = dot(F)
  []
  [cauchy]
    type = CauchyStress
    cauchy_stress = cauchy
    first_piola_kirchhoff_stress = pk1
    deformation_gradient = F
  []
[]

[Postprocessors]
  # Average through-thickness Cauchy stress on each interface (negative = compressive).
  [sig_yy_cu_anode]
    type = SideAverageValue
    variable = sigma_yy
    boundary = 'if_cu_anode'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [sig_yy_anode_sep]
    type = SideAverageValue
    variable = sigma_yy
    boundary = 'if_anode_sep'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [sig_yy_sep_cathode]
    type = SideAverageValue
    variable = sigma_yy
    boundary = 'if_sep_cathode'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [sig_yy_cathode_al]
    type = SideAverageValue
    variable = sigma_yy
    boundary = 'if_cathode_al'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Contact pressure pc = -sigma_yy on each interface [MPa].
  [pc_cu_anode]
    type = ParsedPostprocessor
    pp_names = 'sig_yy_cu_anode'
    expression = '-sig_yy_cu_anode'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [pc_anode_sep]
    type = ParsedPostprocessor
    pp_names = 'sig_yy_anode_sep'
    expression = '-sig_yy_anode_sep'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [pc_sep_cathode]
    type = ParsedPostprocessor
    pp_names = 'sig_yy_sep_cathode'
    expression = '-sig_yy_sep_cathode'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [pc_cathode_al]
    type = ParsedPostprocessor
    pp_names = 'sig_yy_cathode_al'
    expression = '-sig_yy_cathode_al'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Contact resistance per interface [mOhm cm^2].
  [Rc_cu_anode]
    type = ParsedPostprocessor
    pp_names = 'pc_cu_anode'
    expression = '${R0_mOhm_cm2} * exp(-${alpha_per_MPa} * pc_cu_anode)'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Rc_anode_sep]
    type = ParsedPostprocessor
    pp_names = 'pc_anode_sep'
    expression = '${R0_mOhm_cm2} * exp(-${alpha_per_MPa} * pc_anode_sep)'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Rc_sep_cathode]
    type = ParsedPostprocessor
    pp_names = 'pc_sep_cathode'
    expression = '${R0_mOhm_cm2} * exp(-${alpha_per_MPa} * pc_sep_cathode)'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Rc_cathode_al]
    type = ParsedPostprocessor
    pp_names = 'pc_cathode_al'
    expression = '${R0_mOhm_cm2} * exp(-${alpha_per_MPa} * pc_cathode_al)'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Stack-averaged contact pressure and contact-stability metric.
  [pc_avg]
    type = ParsedPostprocessor
    pp_names = 'pc_cu_anode pc_anode_sep pc_sep_cathode pc_cathode_al'
    expression = '0.25 * (pc_cu_anode + pc_anode_sep + pc_sep_cathode + pc_cathode_al)'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Pi_c]
    type = ParsedPostprocessor
    pp_names = 'pc_avg'
    expression = 'pc_avg / ${pmin}'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Temperature trace for the output CSV (constant T_ref during cure, then
  # the DARPA cycle for t > 1).
  [T_pp]
    type = FunctionValuePostprocessor
    function = T_func
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Running min/max of Pi_c over the simulation (useful when soc cycles).
  [Pi_c_min]
    type = TimeExtremeValue
    postprocessor = Pi_c
    value_type = min
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [Pi_c_max]
    type = TimeExtremeValue
    postprocessor = Pi_c
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true
  line_search = none

  l_max_its = 200
  l_tol = 1e-6
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  nl_max_its = 20

  start_time = 0.0
  end_time = ${end_time}
  dt = 0.1
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
  file_base = 'rst/pressure_free_stack'
[]
