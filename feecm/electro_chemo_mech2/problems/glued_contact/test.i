# This is a model with 5 flat blocks
# Block 1 -> Quartz 50 um
# Block 2 -> Cu current collector -> 200 nm = 0.2 um
# Block 3 -> Au bond layer -> 50 nm = 0.05 um
# Block 3 -> Li metal -> 200 um
# Block 4 -> Interface Li layer
# Block 5 -> Ceramic Electrolyte -> LLZ0 -> 500 um (parametrized for curvature measurements)
# Li_ion_conc -- modeled as a steady state thermal problem
# All units in um
# Model is axisymmetric with overall 1000 um Radius -> 1mm (Actual = 1cm)
# This is modified from actual dimensions to simplify simulations
# The contact model is one of glued contact ... where Li does not easily separate from LLZO
# Initial stress applied to Li layer -> Approx 3 MPa to model initial deformation prior to deposition
# Contact resistance is obtained from Juny spread sheet of resistance...
# # Measured ionic-conductivity is assumed to be LLZO = 0.1 mS/cm... which means the overall
# # conductance of a thickness t is K= 0.1 * A/t ---> resistance = 1/K (ohm)
# # Example t = 500 um LLZO, A = pi * 0.5^2 @ 0.1 mS/cm -> K = 0.1e-3 * A/t = 600 ohms
## Measured resistance from simple ohmic calculation (not impedance measurement ) = 10000 ohms
###  Ignoring electronic resistance of metals in the system and the liquid electroylyte
### We can gather that the remaining resistance is coming from the interfaces involved in electrochemistry
 # Possible sources are
 # 1) Interface between Li and liquid electrolyte
 # 2) Interface between Li and LLZO
# In this analysis, since the Liquid interface is completely absent, we make the interface resistance of that interface high
 # Namely voltage drop across that interface is high
 # 1) Rct at liquid-Li interface
 # 2) Rct at Li-LLZO interface = 1000 ohm-cm^2 => R_int = 1000 * pi * 0.5 ^2 ~ 1000 ohms
 # 3) We can also have variation in the Rct at the interface and allow it to vary between 5000 ohm-cm^2 high initial_resistance
 #    and low reistance 1000 ohm-cm^2
 # 4) All other potential drops are assumed to occur at the Li-Liquid interface
 # 5) In order to minimize
# In terms of the units used here -> for very thin film t = 0.01 um, the R =
# --- Mechanical Material properties
# Quartz -> Elastic
#      E = 76 GPa
#      nu = 0.3
# Cu current collector -> Elastic Plastic
#      E = 130 GPa
#      nu = 0.25
#     Yield Strength = 200 MPa
#     Linear Hardening with Hardening modulus = 50 MPa is assumed
# Currently Au bond coat is ignored
# Li- Metal properties are obtained from Wang and Sakamoto
### --- To be filled out
# LLZo - Elastic
#     E = 170 GPa
#     nu = 0.3
# -----------------------------------------------------------------------------------------------------------------------------
xmax = 200
nx = 500
nx2 = 250
[Mesh]
  patch_size = 120
  patch_update_strategy = auto
  displacements = 'ux uy'
  parallel_type = REPLICATED
  [./Quartz_layer]
    type = GeneratedMeshGenerator
    xmax = 500
    xmin = 0.0
    ymax = 70.4
    ymin = 20.4
    nx = 500
    ny = 5
    dim = 2
    elem_type = QUAD4
  [../]
  [./Quartz]
    type = RenameBlockGenerator
    input = 'Quartz_layer'
    old_block_id = '0'
    # new_block_name = 'Quartz'
    new_block_id = '1'
  [../]
  [./Cu_layer]
    type = GeneratedMeshGenerator
    xmax = 500
    xmin = 0.0
    ymax = 20.4
    ymin = 20.2
    nx = 500
    ny = 2
    dim = 2
    elem_type = QUAD4
  [../]
  [./Cu]
    type = RenameBlockGenerator
    input = 'Cu_layer'
    old_block_id = '0'
    new_block_id = '2'
  [../]

  [./Li_layer]
    type = GeneratedMeshGenerator
    xmax = 500
    xmin = 0
    ymax = 20.2
    ymin = 0.2
    nx = 500
    ny = 12
    dim = 2
    bias_y = 1.15
    elem_type = QUAD4
  [../]
  [./Li_1]
    type = RenameBlockGenerator
    input = 'Li_layer'
    old_block_id = '0'
    new_block_id = '3'
  [../]

  [./Li_interphase]
    type = GeneratedMeshGenerator
    xmax = 500
    xmin = 0.0
    ymax = 0.2
    ymin = 0.0
    nx = 500
    ny = 3
    dim = 2
    elem_type = QUAD4
  [../]
  [./Li_2]
    type = RenameBlockGenerator
    input = 'Li_interphase'
    old_block_id = '0'
    new_block_id = '4'
  [../]

  [./upper]
    type = StitchedMeshGenerator
    inputs = 'Quartz Cu Li_1 Li_2'
    stitch_boundaries_pairs = 'bottom top; bottom top; bottom top'
    clear_stitched_boundary_ids = true
  [../]
  [./upper_boundary]
    type = RenameBoundaryGenerator
    input = upper
    old_boundary_name = 'top bottom left right'
    # new_boundary_name = 'upper_top upper_bottom upper_left upper_right'
    new_boundary_id = '10 20 30 40'
  [../]
  [./LLZO_layer_1]
    type = GeneratedMeshGenerator
    dim = 2
    elem_type = QUAD4
    xmax = 500
    xmin = 0.0
    ymax = 0.0
    ymin = -0.5
    nx = 250
    ny = 2
    # boundary_name_prefix = 'LLZO'
  [../]
  [./LLZO_1]
    type = RenameBlockGenerator
    input = 'LLZO_layer_1'
    old_block_id = '0'
    new_block_id = '5'
  [../]
  [./LLZO_layer_2]
    type = GeneratedMeshGenerator
    dim = 2
    elem_type = QUAD4
    xmax = 500
    xmin = 0.0
    ymax = -0.5
    ymin = -50.0
    nx = 250
    ny = 20
    # boundary_name_prefix = 'LLZO'
  [../]
  [./LLZO_2]
    type = RenameBlockGenerator
    input = 'LLZO_layer_2'
    old_block_id = '0'
    new_block_id = '6'
  [../]
  [./LLZO]
    type = StitchedMeshGenerator
    inputs = 'LLZO_1 LLZO_2'
    stitch_boundaries_pairs = 'bottom top'
    clear_stitched_boundary_ids = true
  [../]

  [./LLZO_boundary]
    type = RenameBoundaryGenerator
    input = LLZO
    old_boundary_name = 'top bottom left right'
    # new_boundary_name = 'lower_top lower_bottom lower_left lower_right'
    new_boundary_id = '50 60 70 80'
  [../]

  [./all_mesh]
    type = CombinerGenerator
    inputs = 'upper_boundary LLZO_boundary'
  [../]
  [./secondary_block]
    type = LowerDBlockFromSidesetGenerator
    input = all_mesh
    sidesets = 20
    new_block_id = 200
  [../]

  [./primary_block]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_block
    sidesets = 50
    new_block_id = 500
  [../]
  [./top_center_node]
    type = BoundingBoxNodeSetGenerator
    input = primary_block
    new_boundary = 1000
    location = INSIDE
    bottom_left = '-1 70.3 0'
    top_right = '0.1 70.5 0'
  [../]
  [./bottom_center_node]
    type = BoundingBoxNodeSetGenerator
    input = top_center_node
    new_boundary = 1001
    location = INSIDE
    bottom_left = '-0.1 -50.1 0'
    top_right = '0.1 -49.9 0'
  [../]
  [./left_crack]
    type = BoundingBoxNodeSetGenerator
    input = bottom_center_node
    new_boundary = 1002
    location = INSIDE
    bottom_left = '-0.1 -50.0 0'
    top_right = '0.01 -5.0 0'
  [../]
[]


[GlobalParams]
  displacements = 'ux uy'
[]

[Variables]
  [./ux]
    block = '1 2 3 4 5 6'
  [../]

  [./uy]
    block = '1 2 3 4 5 6'
  [../]

#  [./normal_lm_x]
#    block = 200
#  [../]

#  [./normal_lm_y]
#    block = 200
#  [../]

  # [./thermal_lm]
  #   block = 200
  # [../]
  #
  # [./li_metal_conc]
  #   block = '4'
  #   initial_condition = 0
  # [../]
  #
  # [./li_ion_V]
  #   block = '1 2 3 4 5 6'
  # [../]

[]


[Functions]
  [./pressure]
    type = PiecewiseLinear
    data_file = 'pressure_time1.csv'
    format = columns
  [../]
  [./flux]
    type = ParsedFunction
    # value = 'if (t >= 1, 4e-13, 0.0)'
    value = '4e-4'
  [../]
  [./gapk]
    type = PiecewiseLinear
    data_file = 'gap_cond.csv'
    format = columns
  [../]
  [./gapk1]
    type = PiecewiseLinear
    x = '-10.0 0 1e-2 1 10 100'
    y = '1e-2 1e-2 1e-2 1e-2 0 0'
  [../]

  [./k_function]
    type = ParsedFunction
    # value = 'if (t <= 0.1, 1e-9, 1e-9)'
    value = 'if (t <=0.1, 1e-9, 0.5e-9*cos(2.0*3.14157 * x/10.0) + 0.5e-9 + 1e-11)'
    # value = 'if (t <=0.1, 1e-9, (1.0/(1.0 + exp(2.0*2.5*((x%20.0)-10.0))))*(1e-9-1e-11) + 1e-11)'
    # value = '(1.0/(1.0 + exp(2.0*2.5*(x-20.0))))*(1e-11-1e-13) + 1e-13)'

  [../]
[]

[AuxVariables]
  [./plastic_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_rate]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./li_ion_flux_x]
    order = FIRST
    family = MONOMIAL
    block = '1 2 3 4 5 6'
  [../]
  [./li_ion_flux_y]
    order = FIRST
    family = MONOMIAL
    block = '1 2 3 4 5 6'
  [../]
  [./li_ion_flux_z]
    order = FIRST
    family = MONOMIAL
    block = '1 2 3 4 5 6'
  [../]

  [./bndliflux]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3 4 5 6'
  [../]
[]

[AuxKernels]
  [./bnd_li_flux]
    type = DiffusionFluxNormalToBoundaryAux
    variable = bndliflux
    boundary = '20 50'
    diffusion_variable = li_ion_V
    diffusivity = thermal_conductivity
  [../]
  [./peeq]
    type = ADMaterialRealAux
    variable = plastic_strain
    property = effective_plastic_strain
    execute_on = timestep_end
    block = '3 4'
  [../]
  [./strain_rate]
    type = ADMaterialRealAux
    variable = strain_rate
    property = plastic_strain_rate
    block = '3 4'
  [../]

  [./li_ion_flux_x]
    type = ADDiffusionFluxAux
    variable = li_ion_flux_x
    component = x
    diffusion_variable = li_ion_V
    diffusivity = thermal_conductivity
    block = '1 2 3 4 5 6'
  [../]

  [./li_ion_flux_y]
    type = ADDiffusionFluxAux
    variable = li_ion_flux_y
    component = y
    diffusion_variable = li_ion_V
    diffusivity = thermal_conductivity
    block = '1 2 3 4 5 6'
  [../]
  [./li_ion_flux_z]
    type = ADDiffusionFluxAux
    variable = li_ion_flux_z
    component = z
    diffusion_variable = li_ion_V
    diffusivity = thermal_conductivity
    block = '1 2 3 4 5 6'
  [../]

[]

[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    use_displaced_mesh = true
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress hydrostatic_stress'
    use_automatic_differentiation = true
    block = '1 2 3 4 5 6'
    displacements = 'ux uy'
    extra_vector_tags = 'ref'
  [../]
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'ux uy normal_lm_x normal_lm_y; li_ion_V li_metal_conc thermal_lm'
  acceptable_iterations = 2
[]

[Contact]
  [./mech]
    model = glued
    formulation = kinematic
    penalty = 1e3
    primary = 50
    secondary = 20
  [../]
[]

[Constraints]
  [./thermal_constraint]
    type = GapDisplacementConductanceConstraint
    variable = thermal_lm
    secondary_variable = li_ion_V
    # k_function = gapk1
    k = 1e-3
    use_displaced_mesh = true
    compute_lm_residuals = true
    primary_subdomain = 500
    secondary_subdomain = 200
    primary_boundary = 50
    secondary_boundary = 20
    # extra_vector_tags = 'ref'
  [../]



  # [./normal_lm_x]
  #   type = EqualValueConstraint
  #   variable = normal_lm_x
  #   primary_variable = ux
  #   secondary_variable = ux
  #   primary_subdomain = 500
  #   secondary_subdomain = 200
  #   primary_boundary = 50
  #   secondary_boundary = 20
  #   use_displaced_mesh = true
  # [../]
  #
  # [./normal_lm_y]
  #   type = EqualValueConstraint
  #   variable = normal_lm_y
  #   primary_variable = uy
  #   secondary_variable = uy
  #   primary_subdomain = 500
  #   secondary_subdomain = 200
  #   primary_boundary = 50
  #   secondary_boundary = 20
  #   use_displaced_mesh = true
  # [../]

[]

[Kernels]
  [./li_ion_V]
    type = ADHeatConduction
    block = '1 2 3 4 5 6'
    variable = li_ion_V
    use_displaced_mesh = false
  [../]
  [./li_metal2]
    type = ADChemoMechanoAnsioDiffusion
    block = '4'
    variable = li_metal_conc
    diffusivity = diffusivity
    # diffusivity = diffusivity_Li
    use_displaced_mesh = false
  [../]
  [./li_metal_dt]
    type = ADTimeDerivative
    block = '4'
    variable = li_metal_conc
    use_displaced_mesh = false
  [../]
[]

[Materials]
  # [./diffusivity_Li1]
  #   # type = GenericConstantMaterial
  #   # prop_names = 'diffusivity_Li specific_heat1'
  #   # prop_values = '1e-1 1'
  #   block = '1 2 3 4'
  #   type = ADConstantAnisotropicMobility
  #   tensor = '0 0 0
  #             0 1e15 0
  #             0 0 0'
  #   M_name = 'diffusivity_Li'
  # [../]

  [./diffusivity_Li]
    type = ADDiffusionAlongPrincipalDirectionsMaterial
    diffusivity_vector = '1e-1 0 0'
    block = '4'
  [../]

  [./thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1.0e-2
    block = '6'
  [../]
  [./thermal_function]
    type = ADGenericFunctionMaterial
    prop_names = 'thermal_conductivity'
    prop_values = k_function
    block = '5'
  [../]


  [./thermal_conductivity4]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1e4
    block = '1 2 3 4'
    use_displaced_mesh = true
  [../]

  [./elasticity_tensor_Li]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 7.8e3
    poissons_ratio = 0.3
    block = '3 4'
  [../]
  [./elasticity_tensor_Llzo]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 170e3
    poissons_ratio = 0.25
    block = '5 6'
  [../]
  [./stress_llzo]
    type = ADComputeFiniteStrainElasticStress
    block = '5 6'
  [../]

  [./stress_Li]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    # perform_finite_strain_rotations = true
    block = '4'
  [../]

  [./stress_Li2]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas2'
    # perform_finite_strain_rotations = true
    block = '3'
  [../]
  [./plas]
    type = ADIsoTropicHyperViscoSwellingCreep
    absolute_tolerance = 1e-5
    # relative_tolerance = 1e-06
    hardening_exponent = 2.0
    saturation_resistance = 2.0
    initial_resistance = 0.95
    hardening_modulus = 10.0
    rate_exponent = 0.15
    # reference_strain_rate = 0.05
    activation_energy = 37e3
    gas_constant = 8.314462681 # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = 298

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 0.5
    omega = 2
    alpha = '1 0 0'
    concentration = li_metal_conc
    intBnd = 50
    cref = 0
    block = 4
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    enable = true
  [../]

  [./plas2]
    type = ADIsoTropicHyperViscoCreep
    relative_tolerance = 1e-5
    # absolute_tolerance = 1e-6
    hardening_exponent = 2.0
    saturation_resistance = 2.0
    initial_resistance = 0.95
    hardening_modulus = 10.0
    rate_exponent = 0.15
    # reference_strain_rate = 0.05
    activation_energy = 37e3
    gas_constant = 8.314462681 # kJ/K/mol
    saturation_exponent = 0.05
    pre_factor = 4.25e4
    temperature = 298

    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 0.1
    block = 3
  [../]
  [./elasticity_tensor_quartz]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 76e3
    poissons_ratio = 0.25
    block = '1'
  [../]
  [./stress_quartz]
    type = ADComputeFiniteStrainElasticStress
    block = '1'
  [../]
  [./elasticity_tensor_Cu]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 165e3
    poissons_ratio = 0.3
    block = 2
  [../]
  [./stress_cu]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas_cu'
    block = '2'
  [../]
  [./plas_cu]
    type = ADIsotropicPlasticityStressUpdate
    yield_stress = 130
    hardening_constant = 50
    block = 2
  [../]
[]

[BCs]
  [./Li_top_y]
    type = ADDirichletBC
    boundary = '1001'
    variable = uy
    value = 0.0
    preset = true
  [../]

  # [./Li_top_x]
  #   type = ADDirichletBC
  #   boundary = '10'
  #   variable = ux
  #   value = 0.0
  #   preset = true
  # [../]

  [./left_right]
    type = ADDirichletBC
    preset = true
    variable = ux
    boundary = '30'
    value = 0
  [../]
  [./crack]
    type = ADDirichletBC
    variable = ux
    boundary = 1002
    value = 0
  [../]
  # [./bot_pressure]
  #   type = ADPressure
  #   boundary = '60'
  #   variable = uy
  #   component = 1
  #   function = pressure
  #   use_displaced_mesh = true
  # [../]

  [./li_metal_flux]
    type = ScaledCoupledVarNeumannBC
    boundary = '20'
    variable = li_metal_conc
    v = bndliflux
    scale = -1.036426e-2
    extra_vector_tags = 'ref'
    # value = 1.0
  [../]
  [./li_ion_bot]
    type = FunctionNeumannBC
    boundary = '60'
    variable = li_ion_V
    function = flux
    extra_vector_tags = 'ref'
  [../]
  [./li_ion_top]
    type = ADDirichletBC
    preset = true
    boundary = '10'
    variable = li_ion_V
    value = 0.0
  [../]

[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./bottom_current]
    type = ADSideFluxIntegral
    variable = li_ion_V
    boundary = 50
    diffusivity = thermal_conductivity
  [../]
  [./top_current]
    type = ADSideFluxIntegral
    variable = li_ion_V
    boundary = 10
    diffusivity = thermal_conductivity
  [../]
  [./lagrange]
    type = ElementIntegralVariablePostprocessor
    variable = thermal_lm
    block = 200
  [../]
  # [./auxcurrent]
  #   type = SideIntegralVariablePostprocessor
  #   variable = bndliflux
  #   boundary = 15
  # [../]
  [./over_potential]
    type = SideAverageValue
    variable = li_ion_V
    boundary = 60
  [../]

  [./ext_pressure]
    type = SideAverageValue
    variable = stress_yy
    boundary = 10
  [../]
  [./matl_ts_min]
    type = MaterialTimeStepPostprocessor
    block = '3 4'
  [../]

[]

[Executioner]

  type = Transient

  #Preconditioned JFNK (default)
  solve_type = 'NEWTON'
  automatic_scaling = true
  compute_scaling_once = false
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor'
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1     NONZERO               1e-18               '

  #line_search = contact
  l_max_its = 50
  nl_max_its = 25
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-6
  start_time = 0.0
  dt = 1.0
  dtmax = 5.0
  dtmin = 1e-5
  end_time = 3000.0
  snesmf_reuse_base = true
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.025
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 40
    timestep_limiting_postprocessor = matl_ts_min
  [../]


[] # Executioner

[Outputs]
  # [./check]
  #   type = Checkpoint
  #   num_files = 2
  #   interval = 5
  #   sync_only = true
  #   sync_times = '1 2.0 5.0 10.0 15.0 20.0 25.0 50.0 100.0 200.0'
  # [../]
  [./out]
    type = Exodus
    execute_on =  'INITIAL TIMESTEP_END'
    elemental_as_nodal = true
    # execute_elemental_variables = true
    execute_elemental_on = 'TIMESTEP_END'
    # interval = 10
    sync_only = false
    sync_times = '1.0 15.0 20.0 50.0 80.0 110.0 150.0 180.0 210.0 240.0 270.0 300.0 330.0 360.0 500.0'
    # sync_times = '1 2.0 5.0 10.0 15.0 20.0 25.0 50.0 100.0 200.0'
    # output_material_properties = true
  [../]
  [./checkpoint]
    type = Checkpoint
    num_files = 4
    start_time = 30.0
    sync_times = '30.0 40.0 50.0 70.0 100.0 200.0 250.0 300.0 400.0 500.0'
    sync_only = true
  [../]

  [./csv]
    type = CSV
  [../]
[] # Outputs

[Debug]
  show_var_residual_norms = true
  show_material_props = true
[]
