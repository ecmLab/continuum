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
# ----------------------------------- Note ----------------------------------
# Results in comparison to paper will not be exact since the yield stress is
# not a function of the concentration here. That should be easy to incorporate
# later if necessary
# ----------------------------------------------------------------------------
[Mesh]
  [./mesh]
    type = FileMeshGenerator
    file = 'si_2d_w_se.msh'
  [../]
  [./primary_boundary_anode]
    type = LowerDBlockFromSidesetGenerator
    input = mesh
    sidesets = 'Electrolyte_top'
    new_block_name = 'anode_p_block'
  [../]
  [./secondary_boundary_anode]
    type = LowerDBlockFromSidesetGenerator
    input = primary_boundary_anode
    sidesets = 'Anode_bottom'
    new_block_name = 'anode_s_block'
  [../]
[]

[GlobalParams]
  displacements = 'ux uy'
[]

[Variables]
  [./ux]
    block = 'Electrolyte Anode'
  [../]
  [./uy]
    block = 'Electrolyte Anode'
  [../]
  # -Voltage -
  [./V]
    block = 'Electrolyte Anode'
  [../]

  [./li_metal_conc]
    initial_condition = 614.172e-6
    block = 'Anode'
  [../]
  # -- Stress based chemical potenttial -- #
  [./mu_sigma]
    block = 'Anode'
  [../]
  [./V_anode]
    block = 'anode_s_block'
  [../]
  [./normal_lm_x_anode]
    block = 'anode_s_block'
  [../]
  [./normal_lm_y_anode]
    block = 'anode_s_block'
  [../]

[]


[Constraints]
  [./anode_butler_volmer]
    type = GapDisplacementConductanceConstraint
    variable = V_anode
    secondary_variable = V
    secondary_boundary =  'Anode_bottom'
    secondary_subdomain = 'anode_s_block'
    primary_boundary = 'Electrolyte_top'
    primary_subdomain = 'anode_p_block'
    # k_function = gapk
    k = 1e-2
    use_displaced_mesh = true
    displacements = 'ux uy'
    compute_lm_residuals = true
    include_equilibrium_potential = false
  [../]
  [./anode_mech_x]
    type = EqualValueConstraint
    variable = normal_lm_x_anode
    primary_variable = ux
    secondary_variable = ux
    secondary_boundary =  'Anode_bottom'
    secondary_subdomain = 'secondary_contact_block'
    primary_boundary = 'Electrolyte_top'
    primary_subdomain = 'primary_contact_block'
    use_displaced_mesh = true
  [../]
  [./anode_mech_y]
    type = EqualValueConstraint
    variable = normal_lm_y_anode
    primary_variable = uy
    secondary_variable = uy
    secondary_boundary =  'Anode_bottom'
    secondary_subdomain = 'secondary_contact_block'
    primary_boundary = 'Electrolyte_top'
    primary_subdomain = 'primary_contact_block'
    use_displaced_mesh = true
  [../]
  [./anode_conc]
    type = ScaledBCConstraint
    variable = V_anode
    secondary_variable = li_metal_conc
    primary_variable = V
    primary = false
    secondary_boundary =  'Anode_bottom'
    secondary_subdomain = 'secondary_contact_block'
    primary_boundary = 'Electrolyte_top'
    primary_subdomain = 'primary_contact_block'
    scale = 1.036426e-2
    # extra_vector_tags = 'ref'
  [../]

[]

[AuxVariables]
  [./flux_x]
    order = FIRST
    family = MONOMIAL
    block = 'Electrolyte Anode'
  [../]
  [./flux_y]
    order = FIRST
    family = MONOMIAL
    block = 'Electrolyte Anode'
  [../]
  [./flux_z]
    order = FIRST
    family = MONOMIAL
    block = 'Electrolyte Anode'
  [../]
  [./bndliflux]
    order = CONSTANT
    family = MONOMIAL
    block = 'Electrolyte Anode'
  [../]
  [./growth_potential]
    order = CONSTANT
    family = MONOMIAL
    block = 'Anode'
  [../]
  [./swelling]
    order = CONSTANT
    family = MONOMIAL
    block = 'Anode'
  [../]
  [./elastic_potential]
    order = CONSTANT
    family = MONOMIAL
    block = 'Anode'
  [../]
  [./plastic_strain]
    order = CONSTANT
    family = MONOMIAL
    block = 'Anode'
  [../]
  [./equil_potential]
    order = CONSTANT
    family = MONOMIAL
    block = 'Anode'
  [../]
  [./growth_stress_hydro]
    order = CONSTANT
    family = MONOMIAL
    block = 'Anode'
  [../]
  [./chemical_potential]
    order = CONSTANT
    family = MONOMIAL
    block = 'Anode'
  [../]

[]

[AuxKernels]
  [./growth_stress_hydro]
    type = ADRankTwoScalarAux
    variable = growth_stress_hydro
    rank_two_tensor = growth_stress
    scalar_type = Hydrostatic
    block = 'Anode'

  [../]

  [./chemical_potential]
    type = ADMaterialRealAux
    variable = chemical_potential
    property = chemical_potential
    block = 'Anode'

  [../]

  [./equil_potential]
    type = ADMaterialRealAux
    variable = equil_potential
    property = equilibrium_potential
    block = 'Anode'

  [../]
  [./swelling]
    type = ADMaterialRealAux
    variable = swelling
    property = swelling_vol_change
    block = 'Anode'
  [../]
  [./growth_potential]
    type = ADMaterialRealAux
    variable = growth_potential
    property = swelling_chemical_potential
    block = 'Anode'

  [../]
  [./elastic_potential]
    type = ADMaterialRealAux
    variable = elastic_potential
    property = elastic_chemical_potential
    block = 'Anode'

  [../]
  [./bnd_li_flux]
    type = DiffusionFluxNormalToBoundaryAux
    variable = bndliflux
    boundary = 'Electrolyte_top'
    diffusion_variable = V
    diffusivity = thermal_conductivity
  [../]

  [./li_ion_flux_x]
    type = ADDiffusionFluxAux
    variable = flux_x
    component = x
    diffusion_variable = V
    diffusivity = thermal_conductivity
    block = 'Electrolyte Anode'
  [../]

  [./li_ion_flux_y]
    type = ADDiffusionFluxAux
    variable = flux_y
    component = y
    diffusion_variable = V
    diffusivity = thermal_conductivity
    block = 'Electrolyte Anode'
  [../]
  [./li_ion_flux_z]
    type = ADDiffusionFluxAux
    variable = flux_z
    component = z
    diffusion_variable = V
    diffusivity = thermal_conductivity
    block = 'Electrolyte Anode'
  [../]
  [./plastic_strain]
    type = ADMaterialRealAux
    variable = plastic_strain
    property = effective_plastic_strain
    block = 'Anode'

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
    block = 'Electrolyte Anode'
    # use_finite_deform_jacobian = true
  [../]
[]

[Problem]
  type = ReferenceResidualProblem
  # solution_variables = 'ux uy normal_lm thermal_lm li_ion_V thermal_lm2 li_metal_conc'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'ux uy; li_metal_conc mu_sigma'
  acceptable_iterations = 2
  # restart_file_base = check/full_model_cp/0002
[]

[ICs]
  [./v_anode]
    type = ConstantIC
    value = 1000
    variable = V
    block = 'Anode'
  [../]
[]

[Kernels]
  [./li_ion_V]
    type = ADHeatConduction
    variable = V
    use_displaced_mesh = false
    block = 'Electrolyte Anode'
  [../]
  [./stress_chemical_potential]
    type = ADStressBasedChemicalPotential
    variable = mu_sigma
    use_displaced_mesh = false
    extra_vector_tags = 'ref'
    block = 'Anode'
  [../]
  [./Diffusion]
    type = ADChemoMechanoAnsioDiffusion
    variable = li_metal_conc
    # stress_based_chemical_potential = mu_sigma
    diffusivity = diffusivity
    block = 'Anode'
    # R = 96.4853329
    # temperature = 298
  [../]
  [./diffusion_dt]
    type = ADTimeDerivative
    variable = li_metal_conc
    block = 'Anode'
  [../]

[]

[Materials]
  # [./diffusivity_Li]
  #   type = ADIsotropicDiffusionMaterial
  #   diffusion_coef = 3e-7
  # [../]
  [./diffusivity_Li]
    type = ADDiffusionAlongPrincipalDirectionsMaterial
    diffusivity_vector = '3e-3 0 0'
    block = 'Anode'
  [../]

  [./equilibrium_potential]
    type = ADComputeEquilibriumPotential
    R = 8.31446
    faraday = 96.4853329
    temperature = 298
    cref = 7.874e-2
    concentration = li_metal_conc
    include_conc = true
    include_reaction_rate = true
    # reaction_rate = 780.0
    reaction_rate_function = reaction_rate
    include_mechanical_effects = true
    exclude_elastic_contribution = false
    block = 'Anode'
  [../]


  [./thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1
    block = 'Anode'
  [../]
  [./thermal_conductivity2]
    type = ADHeatConductionMaterial
    thermal_conductivity = 0.1
    block = 'Electrolyte'
  [../]

  [./elasticity_tensor_Si]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 100.0e3
    poissons_ratio = 0.26
    block = 'Anode Electrolyte'
  [../]
  [./stress_llzo]
    type = ADComputeFiniteStrainElasticStress
    block = 'Electrolyte'
  [../]
  [./stress_Si]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    block = 'Anode'
  [../]
  [./plas]
    type = ADIsoTropicHyperViscoSwelling
    absolute_tolerance = 1e-5
    relative_tolerance = 1e-8
    isotropic_swelling = false
    rate_form = false
    # alpha = '0.33333333 0.33333333 0.33333333'
    alpha = '1 0 0'

    hardening_exponent = 1.0
    saturation_resistance = 400.0
    initial_resistance = 120.0
    hardening_modulus = 100.0
    rate_exponent = 0.25
    reference_strain_rate = 0.6e-9
    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 1.0
    omega = 8.889
    concentration = li_metal_conc
    cref = 614.176e-6
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    intBnd = 'Electrolyte_top'
    block = 'Anode'
  [../]
[]
[Functions]
  # --- This is the OCV ---
  # As in Bower (2014) this is a linear function

  [./reaction_rate]
    type = PiecewiseLinear
    x = '0.0 4.857'
    y = '780.0 0.0'
  [../]
[]


[BCs]
  # [./current]
  #   type = ADButlerVolmerBC
  #   current_density = 1.2e-3
  #   exchange_current_density = 1e-6
  #   faraday = 96.4853329
  #   Temperature = 298
  #   R = 8.314462681
  #   variable = V
  #   boundary = top
  #   extra_vector_tags = 'ref'
  # [../]
  [./current]
    type = ADNeumannBC
    variable = V
    boundary = 'Electrolyte_bottom'
    value = 1e-4
    extra_vector_tags = 'ref'
  [../]
  [./disp_x]
    type = ADDirichletBC
    preset = true
    value = 0
    variable = ux
    boundary = 'Electrolyte_bottom'
  [../]
  [./disp_y]
    type = ADDirichletBC
    preset = true
    value = 0
    variable = uy
    boundary = 'Electrolyte_bottom'
  [../]
  [./pressure]
    type = ADPressure
    variable = uy
    component = 1
    function = '1e-6'
    boundary = 'Anode_top'
  [../]
  # [./conc]
  #   type = ADNeumannBC
  #   variable = li_metal_conc
  #   # v = bndliflux
  #   value = 1.24371375e-5
  #   boundary = top
  #   extra_vector_tags = 'ref'
  # [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./Voltage]
    type = SideAverageValue
    boundary = 'Anode_top'
    variable = V
  [../]
  [./stress]
    type = ADElementAverageMaterialProperty
    mat_prop = stress_xx
    block = 'Anode'
  [../]
  [./V0]
    type = SideAverageValue
    boundary = 'Anode_top'
    variable = equil_potential
  [../]
  [./chem_potential]
    type = SideAverageValue
    boundary = 'Anode_top'
    variable = chemical_potential
  [../]
  [./conc]
    type = SideAverageValue
    boundary = 'Anode_top'
    variable = li_metal_conc
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  automatic_scaling = true
  compute_scaling_once = true
  # resid_vs_jac_scaling_param = 0.5
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration'# -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1'#     NONZERO               1e-20               '
  dt = 200
  # num_steps = 20
  nl_max_its = 35
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-6
  dtmax = 200
  verbose = true
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 50
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 40
    # timestep_limiting_postprocessor = matl_ts_min
  [../]
  end_time = 12000
  # num_steps = 20
  # snesmf_reuse_base = true
  # scaling_group_variables = 'ux uy; V li_metal_conc'

[]

[Outputs]
  exodus = true
  csv = true
  [./check]
    type = Checkpoint
    num_files = 2
    interval = 10
  [../]
[]
