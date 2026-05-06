# Test for Si to compare to Bower and Guduru (2014)
# Material Properties
# Molar density of Si = c_0 = 7.874e4 mol/m^3 = 7.874e-2
# Modulus of Si = 100 GPa = 100e3
# poissons_ratio of Si = 0.26
# Vol. expansion of Si (normalized conc ) = 0.7
#    -> Vol. expansion of Si -> 0.7 * 7.874e-2 = omega
# epsilon_0 (strain rate Si) = 0.6e-9 /s
# initial yield stress of Si = 0.12 GPa = 0.12e3
# Stress exponent for plastic flow in Si = 4 -> 0.25
# initial conc of Li in Si = 0.0078 = 0.0078 * 7.87e-2
#applied current density = 0.012 mA/cm^2 = 1.2e-5
# exchange_current_density = 0.001 mA/cm^2 = 1e-6
[Mesh]
  [./Si]
    type = GeneratedMeshGenerator
    nx = 1
    ny = 50
    xmin = 0
    xmax = 0.01
    ymin = 0
    ymax = 0.2
    elem_type = QUAD4
    dim = 2
  [../]
[]
[GlobalParams]
  displacements = 'ux uy'
[]

[Variables]
  [./ux]
  [../]
  [./uy]
  [../]
  [./V]
  [../]
  [./li_metal_conc]
    initial_condition = 614.172e-6
  [../]

[]

[AuxVariables]
  [./V0]
  [../]
  [./flux_x]
    order = FIRST
    family = MONOMIAL
  [../]
  [./flux_y]
    order = FIRST
    family = MONOMIAL
  [../]
  [./flux_z]
    order = FIRST
    family = MONOMIAL
  [../]
  [./bndliflux]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  # [./li_metal_conc]
  #   type = ConstantAux
  #   value = 0
  #   variable = li_metal_conc
  # [../]
  [./bnd_li_flux]
    type = DiffusionFluxNormalToBoundaryAux
    variable = bndliflux
    boundary = 'top'
    diffusion_variable = V
    diffusivity = thermal_conductivity
  [../]

  [./V0]
    type = ConstantAux
    value = 780.0
    variable = V0
  [../]
  [./li_ion_flux_x]
    type = ADDiffusionFluxAux
    variable = flux_x
    component = x
    diffusion_variable = V
    diffusivity = thermal_conductivity
  [../]

  [./li_ion_flux_y]
    type = ADDiffusionFluxAux
    variable = flux_y
    component = y
    diffusion_variable = V
    diffusivity = thermal_conductivity
  [../]
  [./li_ion_flux_z]
    type = ADDiffusionFluxAux
    variable = flux_z
    component = z
    diffusion_variable = V
    diffusivity = thermal_conductivity
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress hydrostatic_stress'
    use_automatic_differentiation = true
    extra_vector_tags = 'ref'
    # use_finite_deform_jacobian = true
  [../]
[]

[Problem]
  type = ReferenceResidualProblem
  # solution_variables = 'ux uy normal_lm thermal_lm li_ion_V thermal_lm2 li_metal_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'ux uy; V li_metal_conc'
  acceptable_iterations = 2
  # restart_file_base = check/full_model_cp/0002
[]


[Kernels]
  [./li_ion_V]
    type = ADHeatConduction
    variable = V
    use_displaced_mesh = false
  [../]
  [./li_metal2]
    type = ADMatDiffusion
    variable = li_metal_conc
    diffusivity = diffusivity
    use_displaced_mesh = false
  [../]
  [./li_metal_dt]
    type = ADTimeDerivative
    variable = li_metal_conc
    use_displaced_mesh = false
  [../]

[]

[Materials]
  [./diffusivity_Li]
    type = ADGenericConstantMaterial
    prop_names = 'diffusivity'
    prop_values = '3e-7'
  [../]
  [./thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1
  [../]
  [./elasticity_tensor_Si]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 100.0e3
    poissons_ratio = 0.26
  [../]
  [./stress_Si]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
  [../]
  [./plas]
    type = ADIsoTropicHyperViscoSwelling
    absolute_tolerance = 1e-8
    relative_tolerance = 1e-8
    isotropic_swelling = true
    rate_form = false
    alpha = '0.33333333 0.33333333 0.33333333'

    hardening_exponent = 1.0
    saturation_resistance = 200.0
    initial_resistance = 120.0
    hardening_modulus = 100.0
    rate_exponent = 0.25
    reference_strain_rate = 0.6e-9
    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 0.1
    omega = 8.889
    concentration = li_metal_conc
    cref = 614.176e-6
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
  [../]
[]

[BCs]
  [./current]
    type = ADButlerVolmerBC
    current_density = 1.2e-5
    exchange_current_density = 1e-6
    faraday = 96.4853
    Temperature = 298
    R = 8.314462681
    equilibrium_potential = V0
    variable = V
    boundary = top
    extra_vector_tags = 'ref'
  [../]
  [./disp_x]
    type = ADDirichletBC
    preset = true
    value = 0
    variable = ux
    boundary = 'left right bottom'
  [../]
  [./disp_y]
    type = ADDirichletBC
    preset = true
    value = 0
    variable = uy
    boundary = 'bottom'
  [../]
  [./pressure]
    type = ADPressure
    variable = uy
    component = 1
    function = '1'
    boundary = 'top'
  [../]
  [./conc]
    type = ScaledCoupledVarNeumannBC
    variable = li_metal_conc
    v = bndliflux
    scale = -1.036428e-2
    boundary = top
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
  compute_scaling_once = false
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration'# -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1'#     NONZERO               1e-20               '
  dt = 200
  # num_steps = 20
  nl_max_its = 35
  nl_abs_tol = 6e-13
  nl_rel_tol = 1e-6
  dtmax = 50
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 50
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 40
    # timestep_limiting_postprocessor = matl_ts_min
  [../]
  end_time = 8000
  snesmf_reuse_base = false
  # scaling_group_variables = 'ux uy; V li_metal_conc'

[]

[Outputs]
  exodus = true
[]
