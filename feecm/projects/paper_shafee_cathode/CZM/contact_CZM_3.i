## --- Bulk Material Properties ---
ymod=50000.00           # Young Modulus of SE [MPa]
pr_se=0.26              # Poissons Ratio of SE
ystr=20                 # Yield Strength of SE [MPa]
ustr=24                 # Ultimate Strength of SE [MPa]
plstr=0.004             # Plastic Strain of SE

ustr_am=12600.0             # Ultimate Strength of AM [MPa]
ymod_am=177500          # Young Modulus of AM [MPa]
pr_am=0.26              # Poissons Ratio of A

sptop=5                 # Stack Pressure [MPa]
alpha_nmc = 0.0002922   # Thermal expansion coefficient of NMC

## --- CZM Parametric Variables ---
cohesion_multiplier = 1.0       # 1.0 = Baseline, 0.5 = 50% Cohesion, 0.0 = No Cohesion
interface_thickness = 0.01      # Effective interface layer thickness [mm] tied to target element edge length (H_EL)

czm_normal_strength_base = 10.0      # Base Normal strength [MPa]
czm_GIc_base             = 0.01      # Base Mode I fracture energy [N/mm]

czm_penalty = ${fparse ymod / interface_thickness}

czm_normal_strength = ${fparse czm_normal_strength_base * cohesion_multiplier}
czm_GIc             = ${fparse czm_GIc_base * cohesion_multiplier}

czm_shear_strength  = ${fparse czm_normal_strength / sqrt(3)}
czm_GIIc            = ${fparse czm_GIc * (czm_shear_strength / czm_normal_strength)^2}


[Problem]
  type = FEProblem
  solve = true
[]

[Mesh]
  [./fmg]
    type = FileMeshGenerator
    # Mesh file in mm
    file = mesh_czm.msh
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
  [./hf_NMC]
    type = PiecewiseLinear
    x = '0'
    y = '${ustr_am}'
  [../]
  [./hf_LPS]
    type = PiecewiseLinear
    x = '0 ${plstr}'
    # Hardening Function LPS
    y = '${ystr} ${ustr}'
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
    youngs_modulus = ${ymod}
    poissons_ratio = ${pr_se}
  [../]
  [./isotropic_plasticity_LPS]
    type = IsotropicPlasticityStressUpdate
    block = 'block_LPS'
    yield_stress = ${ystr}
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

  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu       NONZERO               1e-10'

  line_search = none
  nl_max_its  = 99
  nl_rel_tol  = 1e-6
  nl_abs_tol  = 1e-8
  l_tol       = 1e-8

  start_time = 0.0
  end_time   = 1
  dtmin      = 1e-6

  [TimeStepper]
    type             = IterationAdaptiveDT
    dt               = 0.02
    optimal_iterations = 50
    iteration_window   = 2
    growth_factor    = 1.2
    cutback_factor   = 0.5
  []
[]

[Outputs]
  exodus = true
  file_base = rst2/ystr${ystr}_E${ymod}_spTop${sptop}_CZM_fig6
  [./csv]
    type = CSV
    execute_on = 'final'
  [../]
[]

[Postprocessors]
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
    boundary = 'block_LPS_block_NMC'
    sort_by = x
    execute_on = 'final'
  [../]
[]