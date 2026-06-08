## --- Bulk Material Properties ---
ymod_se=1000             # Young Modulus of SE [MPa] (Sweeping from 100 to 1000 MPa)
Hv_se=30.00              # Vickers Hardness of SE [MPa] (Sweeping from 3 to 30 MPa)
pr_se=0.30              # Poissons Ratio of SE

ustr_am=10000           # Ultimate Strength of AM [MPa] (Some High Value to Ignore AM Plasticity)
ymod_am=90000           # Young Modulus of AM [MPa]
pr_am=0.26              # Poissons Ratio of AM

## --- Calculated Material Properties ---
ustr_se=${fparse Hv_se / 3}                             # Ultimate Strength of SE [MPa] (Coe. = 3)
ystr_se=${fparse ustr_se / 1.2}                         # Yield Strength of SE [MPa]
plstr=${fparse (ustr_se - ystr_se) / (ymod_se / 10)}    # Plastic Strain of SE

## --- Boundary Conditions Properties ---
sptop=10                      # Stack Pressure [MPa]
alpha_nmc=2.9927418e-4        # Thermal expansion coefficient of NMC (8.3% expansion)

## --- CZM Parametric Variables ---
czm_B = 1       # 1.0 = Baseline, 0.0 = No Cohesion (Sweeping from 1 to 10)

interface_thickness = 0.02      # Effective interface layer thickness [mm] tied to target element edge length (H_IFACE)

czm_CED = 10.0           # Cohesion Energy Density [MPa]
czm_GIc_base = 0.01      # Base Mode I fracture energy [N/mm] (0.01 N/mm = 10 J/m^2)

czm_penalty = ${fparse ymod_se / interface_thickness}

czm_normal_strength = ${fparse czm_CED * czm_B}
czm_GIc             = ${fparse czm_GIc_base * czm_B}

czm_shear_strength  = ${fparse czm_normal_strength / sqrt(3)}
czm_GIIc            = ${fparse czm_GIc * (czm_shear_strength / czm_normal_strength)^2}

## --- Output Control ---
output_times = '0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20
                0.22 0.24 0.26 0.28 0.30 0.32 0.34 0.36 0.38 0.40
                0.42 0.44 0.46 0.48 0.50 0.52 0.54 0.56 0.58 0.60
                0.62 0.64 0.66 0.68 0.70 0.72 0.74 0.76 0.78 0.80
                0.82 0.84 0.86 0.88 0.90 0.92 0.94 0.96 0.98 1.00'


[Problem]
  type = FEProblem
  solve = true
[]

[Mesh]
  [./fmg]
    type = FileMeshGenerator
    # Mesh file in mm
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
    block = 'block_NMC'
  [../]
  [./eigenstrain_yy]
    order = FIRST
    family = MONOMIAL
    block = 'block_NMC'
  [../]
  [./total_strain_xx]
    order = FIRST
    family = MONOMIAL
    block = 'block_NMC'
  [../]
  [./total_strain_yy]
    order = FIRST
    family = MONOMIAL
    block = 'block_NMC'
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
  [./hf_NMC]
    type = PiecewiseLinear
    x = '0'
    y = '${ustr_am}'
  [../]
  [./hf_LPS]
    type = PiecewiseLinear
    x = '0 ${plstr}'
    # Hardening Function LPS
    y = '${ystr_se} ${ustr_se}'
  [../]
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [./NMC]
        strain = FINITE
        add_variables = true
        eigenstrain_names = eigenstrain
        generate_output = 'stress_xx stress_yy vonmises_stress strain_xx strain_yy'
        block = 'block_NMC'
      [../]
      [./LPS]
        strain = FINITE
        add_variables = true
        generate_output = 'stress_yy stress_xx vonmises_stress plastic_strain_xx plastic_strain_yy'
        block = 'block_LPS'
      [../]
    [../]
  [../]
[]

[Physics/SolidMechanics/CohesiveZone]
  [./czm_arc]
    strain = FINITE
    boundary = 'block_NMC_block_LPS'
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
    type = RankTwoAux
    block = 'block_NMC'
    rank_two_tensor = eigenstrain
    variable = eigenstrain_yy
    index_i = 1
    index_j = 1
    execute_on = 'initial timestep_end'
  [../]
  [./eigenstrain_xx]
    type = RankTwoAux
    block = 'block_NMC'
    rank_two_tensor = eigenstrain
    variable = eigenstrain_xx
    index_i = 0
    index_j = 0
    execute_on = 'initial timestep_end'
  [../]
  [./total_strain_yy]
    type = RankTwoAux
    block = 'block_NMC'
    rank_two_tensor = total_strain
    variable = total_strain_yy
    index_i = 1
    index_j = 1
    execute_on = 'initial timestep_end'
  [../]
  [./total_strain_xx]
    type = RankTwoAux
    block = 'block_NMC'
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
  [./elasticity_tensor_NMC]
    type = ComputeIsotropicElasticityTensor
    block = 'block_NMC'
    youngs_modulus = ${ymod_am}
    poissons_ratio = ${pr_am}
  [../]
  [./isotropic_plasticity_NMC]
    type = IsotropicPlasticityStressUpdate
    block = 'block_NMC'
    yield_stress = ${ustr_am}
    hardening_function = hf_NMC
  [../]
  [./radial_return_stress_NMC]
    type = ComputeMultipleInelasticStress
    tangent_operator = nonlinear
    inelastic_models = 'isotropic_plasticity_NMC'
    block = 'block_NMC'
  [../]
  [./thermal_expansion_strain_NMC]
    type = ComputeThermalExpansionEigenstrain
    block = 'block_NMC'
    stress_free_temperature = 0
    thermal_expansion_coeff = ${alpha_nmc}
    temperature = temp
    eigenstrain_name = eigenstrain
  [../]

  [./elasticity_tensor_LPS]
    type = ComputeIsotropicElasticityTensor
    block = 'block_LPS'
    youngs_modulus = ${ymod_se}
    poissons_ratio = ${pr_se}
  [../]
  [./isotropic_plasticity_LPS]
    type = IsotropicPlasticityStressUpdate
    block = 'block_LPS'
    yield_stress = ${ystr_se}
    hardening_function = hf_LPS
  [../]
  [./radial_return_stress_LPS]
    type = ComputeMultipleInelasticStress
    tangent_operator = nonlinear
    inelastic_models = 'isotropic_plasticity_LPS'
    block = 'block_LPS'
  [../]

  [./czm_damage]
    type = BiLinearMixedModeTraction
    boundary = 'block_NMC_block_LPS'

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
  nl_abs_tol = 1e-8
  l_tol = 1e-8
  start_time = 0.0
  n_startup_steps = 1
  end_time = 1
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    optimal_iterations = 50
    iteration_window  = 10
    growth_factor = 2.0
    cutback_factor = 0.7
    cutback_factor_at_failure = 0.5
    linear_iteration_ratio = 100
  []
[]

[Outputs]
  file_base = rst_czm/E${ymod_se}_H${Hv_se}_C${czm_B}_spTop${sptop}
  [exodus]
    type = Exodus
    sync_times = '${output_times}'
    sync_only = true
  []
  [csv]
    type = CSV
    sync_times = '${output_times}'
    sync_only = true
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
    boundary = 'block_NMC_block_LPS'
    sort_by = x
    execute_on = 'TIMESTEP_END'
    outputs = csv
  [../]
[]