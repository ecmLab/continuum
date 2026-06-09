## ============================================================================
## UNIT SYSTEM: um (length) - MPa (stress) - uN (force).  1 MPa = 1 uN/um^2.
##   Representative cell = 5 x 5 um; NVP particle radius ~4 um.  The mesh is
##   unit-agnostic (coords read 5.0); these material params set the um scale.
##   Energy release rate G_Ic in MPa*um (= uN/um).
## Catholyte / solid electrolyte = NACS (Na2Al1.35Cl4.05S).  Blocks & vars are
##   named *_NACS / *_nacs (formerly LPS/_se -- a different material).
## ============================================================================
## --- Catholyte (SE) elasto-viscoplastic properties ---
## Rate-dependent (viscoplastic) catholyte: linear elastic + static J2 yield
## (yield from Vickers hardness, Tabor sigma_y ~ H_v/3) + Norton power-law
## creep representing the viscous flow of the plasticized glassy NACS.
##
## ASSUMED temperature presets -- PLACEHOLDERS to be RECALIBRATED to the redone
## hold-segment nanoindentation (25/60/90 C). creep_coef A is tied to the
## (currently proxy) loading timescale and MUST be refit once the real
## charge/discharge loading rate is set.
##         |  25 C rigid | 60 C soft | 90 C plastic/slurry
##   E MPa |    1000     |    500    |    100
##   H_v   |     30      |     12    |     3
##   A     |    1e-7     |    1e-5   |    1e-3
##   n     |     4       |     4     |     4
## Active set = 60 C (override on the command line to sweep T / build the map).
ymod_nacs=500              # Young Modulus of SE [MPa]
Hv_nacs=12.00              # Vickers Hardness of SE [MPa]
creep_coef=1e-5          # Norton creep leading coefficient A [MPa^-n / time] (ASSUMED)
creep_exp=4              # Norton creep stress exponent n (ASSUMED; cf. LPSCl ref n=4)
pr_nacs=0.30              # Poissons Ratio of SE

ymod_nvp=90000           # Young Modulus of NVP [MPa]
pr_nvp=0.26              # Poissons Ratio of NVP

## --- Calculated Material Properties ---
ustr_nacs=${fparse Hv_nacs / 3}                             # Ultimate Strength of SE [MPa] (Coe. = 3)
ystr_nacs=${fparse ustr_nacs / 1.2}                         # Yield Strength of SE [MPa]
plstr=${fparse (ustr_nacs - ystr_nacs) / (ymod_nacs / 10)}    # Plastic Strain of SE

## --- Boundary Conditions Properties ---
sptop=10                      # Stack Pressure [MPa]
alpha_nvp=2.9927418e-4        # Thermal expansion coefficient of NVP (8.3% expansion)

## --- CZM Parametric Variables ---
czm_B = 1       # 1.0 = Baseline, 0.0 = No Cohesion (Sweeping from 1 to 10)

interface_thickness = 0.02      # Effective interface layer thickness [um] tied to target element edge length (H_IFACE)

czm_CED = 10.0           # Cohesion Energy Density [MPa]
czm_GIc_base = 10.0       # Base Mode I fracture energy [MPa*um] (10 MPa*um = 10 J/m^2)

czm_penalty = ${fparse ymod_nacs / interface_thickness}

czm_normal_strength = ${fparse czm_CED * czm_B}
czm_GIc             = ${fparse czm_GIc_base * czm_B}

czm_shear_strength  = ${fparse czm_normal_strength / sqrt(3)}
czm_GIIc            = ${fparse czm_GIc * (czm_shear_strength / czm_normal_strength)^2}

## --- Output Control ---
output_times = '0.05 0.5 1.0'


[Problem]
  type = FEProblem
  solve = true
[]

[Mesh]
  [./fmg]
    type = FileMeshGenerator
    # Mesh file (coords unit-agnostic; um system)
    file = mesh_czm_5050.msh
  [../]
  [./split]
    type = BreakMeshByBlockGenerator
    input = fmg
    split_interface = true
  [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[AuxVariables]
  [./temp]
    initial_condition = 0
  [../]
  [./eigenstrain_xx]
    order = FIRST
    family = MONOMIAL
    block = 'block_NVP'
  [../]
  [./eigenstrain_yy]
    order = FIRST
    family = MONOMIAL
    block = 'block_NVP'
  [../]
  [./total_strain_xx]
    order = FIRST
    family = MONOMIAL
    block = 'block_NVP'
  [../]
  [./total_strain_yy]
    order = FIRST
    family = MONOMIAL
    block = 'block_NVP'
  [../]
[]

[Functions]
  [./temperature_load]
    type = ParsedFunction
    expression = 'if(t<=0.5, 180*t, (-1)*(t-0.5)*180 + 90.0)'
  [../]
  [./press_ramp]
    type = PiecewiseLinear
    x = '0    0.05  1.0'
    y = '0.0  1.0   1.0'
  [../]
  [./hf_NACS]
    type = PiecewiseLinear
    x = '0 ${plstr}'
    # Hardening Function NACS
    y = '${ystr_nacs} ${ustr_nacs}'
  [../]
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [./NVP]
        strain = FINITE
        add_variables = true
        eigenstrain_names = eigenstrain
        use_automatic_differentiation = true
        generate_output = 'stress_xx stress_yy vonmises_stress strain_xx strain_yy'
        block = 'block_NVP'
      [../]
      [./NACS]
        strain = FINITE
        add_variables = true
        use_automatic_differentiation = true
        generate_output = 'stress_yy stress_xx vonmises_stress plastic_strain_xx plastic_strain_yy'
        block = 'block_NACS'
      [../]
    [../]
  [../]
[]

[Physics/SolidMechanics/CohesiveZone]
  [./czm_arc]
    strain = FINITE
    boundary = 'block_NVP_block_NACS'
    generate_output = 'traction_x traction_y normal_traction tangent_traction jump_x jump_y normal_jump tangent_jump'
  [../]
[]

[AuxKernels]
  [./tempfuncaux]
    type = FunctionAux
    variable = temp
    function = temperature_load
  [../]
  [./eigenstrain_yy]
    type = ADRankTwoAux
    block = 'block_NVP'
    rank_two_tensor = eigenstrain
    variable = eigenstrain_yy
    index_i = 1
    index_j = 1
    execute_on = 'initial timestep_end'
  [../]
  [./eigenstrain_xx]
    type = ADRankTwoAux
    block = 'block_NVP'
    rank_two_tensor = eigenstrain
    variable = eigenstrain_xx
    index_i = 0
    index_j = 0
    execute_on = 'initial timestep_end'
  [../]
  [./total_strain_yy]
    type = ADRankTwoAux
    block = 'block_NVP'
    rank_two_tensor = total_strain
    variable = total_strain_yy
    index_i = 1
    index_j = 1
    execute_on = 'initial timestep_end'
  [../]
  [./total_strain_xx]
    type = ADRankTwoAux
    block = 'block_NVP'
    rank_two_tensor = total_strain
    variable = total_strain_xx
    index_i = 0
    index_j = 0
    execute_on = 'initial timestep_end'
  [../]
[]

[BCs]
  [./x_disp]
    type = DirichletBC
    variable = disp_x
    boundary = 'block_left block_right'
    value = 0.0
  [../]
  [./y_disp]
    type = DirichletBC
    variable = disp_y
    boundary = 'block_bottom'
    value = 0.0
  [../]
  [./top_press]
    type = Pressure
    variable = disp_y
    boundary = 'block_top'
    factor   = ${sptop}
    function = press_ramp
  [../]
[]

[Materials]
  [./elasticity_tensor_NVP]
    type = ADComputeIsotropicElasticityTensor
    block = 'block_NVP'
    youngs_modulus = ${ymod_nvp}
    poissons_ratio = ${pr_nvp}
  [../]
  [./elastic_stress_NVP]
    # NVP cathode = isotropic linear elastic (no plasticity).
    type = ADComputeFiniteStrainElasticStress
    block = 'block_NVP'
  [../]
  [./thermal_expansion_strain_NVP]
    type = ADComputeThermalExpansionEigenstrain
    block = 'block_NVP'
    stress_free_temperature = 0
    thermal_expansion_coeff = ${alpha_nvp}
    temperature = temp
    eigenstrain_name = eigenstrain
  [../]

  [./elasticity_tensor_NACS]
    type = ADComputeIsotropicElasticityTensor
    block = 'block_NACS'
    youngs_modulus = ${ymod_nacs}
    poissons_ratio = ${pr_nacs}
  [../]
  [./isotropic_plasticity_NACS]
    type = ADIsotropicPlasticityStressUpdate
    block = 'block_NACS'
    yield_stress = ${ystr_nacs}
    hardening_function = hf_NACS
  [../]
  [./power_law_creep_NACS]
    # Norton viscous flow of the plasticized catholyte (rate dependence).
    # creep_rate = coefficient * vonMises^n   (activation_energy=0 -> Arrhenius
    # term = 1; the T-dependence is baked into coefficient via the preset above).
    type = ADPowerLawCreepStressUpdate
    block = 'block_NACS'
    coefficient = ${creep_coef}        # Norton A (ASSUMED; recalibrate)
    n_exponent  = ${creep_exp}         # Norton n (ASSUMED)
    activation_energy = 0              # no temperature coupling
    max_inelastic_increment = 1e-3     # bound creep increment -> auto substepping
  [../]
  [./radial_return_stress_NACS]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'isotropic_plasticity_NACS power_law_creep_NACS'
    block = 'block_NACS'
  [../]

  [./czm_damage]
    type = BiLinearMixedModeTraction
    boundary = 'block_NVP_block_NACS'

    normal_strength = ${czm_normal_strength}
    shear_strength  = ${czm_shear_strength}
    GI_c  = ${czm_GIc}
    GII_c = ${czm_GIIc}
    penalty_stiffness = ${czm_penalty}

    mixed_mode_criterion = POWER_LAW
    eta = 2.0
    viscosity = 1e-4
  [../]
[]

[Executioner]
  type = Transient
  automatic_scaling = true
  solve_type = NEWTON
  dtmin = 1e-6
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type -mat_mumps_icntl_24'
  petsc_options_value = 'lu       mumps                       1'
  line_search = none
  nl_max_its = 99
  nl_rel_tol = 1e-6
  # abs_tol loosened 1e-8 -> 1e-7: the force-balance residual is well converged
  # by ~1e-7 (Newton then crawls linearly because the CZM/creep tangent is
  # inexact), so 1e-8 only buys wasted iterations. ~2.4x fewer NL iters; gap and
  # traction unchanged to ~1e-5 relative. Use 1e-6 for ~4x if a sweep tolerates
  # slightly looser convergence.
  nl_abs_tol = 1e-7
  l_tol = 1e-8
  start_time = 0.0
  n_startup_steps = 1
  end_time = 1
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    optimal_iterations = 8
    iteration_window  = 2
    growth_factor = 2.0
    cutback_factor = 0.7
    cutback_factor_at_failure = 0.5
    linear_iteration_ratio = 100
  []
[]

[Outputs]
  file_base = rst_czm/E${ymod_nacs}_H${Hv_nacs}_C${czm_B}_spTop${sptop}
  [exodus]
    type = Exodus
    sync_times = '${output_times}'
  []
  [csv]
    type = CSV
    sync_times = '${output_times}'
  []
[]
[Postprocessors]
  [./Gap]
    type = PointValue
    point = '0.0 3.98 0.0'
    variable = disp_y
  [../]
  [./Temp]
     type = ElementAverageValue
     variable = temp
   [../]
[]
[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]
[VectorPostprocessors]
  [./disp_xy_along_arc]
    type = NodalValueSampler
    variable = 'disp_x disp_y'
    # Updated boundary name to the inner sideset generated by BreakMeshByBlockGenerator
    boundary = 'block_NVP_block_NACS'
    sort_by = x
    execute_on = 'TIMESTEP_END'
    outputs = csv
  [../]
[]