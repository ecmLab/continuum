# Test for AgLi problem
# maximum conc of Li in Ag (AgLi9) = 0.0528 mol/cm^3
# initial conc of Li in Ag = 0.0078 (0.0078 * 0.0528 mol/cm^3)
#applied current density = 0.012 mA/cm^2 = 1.2e-5
# exchange_current_density = 0.001 mA/cm^2 = 1e-6
# ----------------------------------- Note ----------------------------------
[Mesh]
  [./Ag]
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

[Variables]
  # -Voltage -
  [./V]
  [../]

  [./li_metal_conc]
    initial_condition = 4.1184e-3
  [../]
[]

[AuxVariables]
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
  [./equil_potential]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./chemical_potential]
    order = CONSTANT
    family = MONOMIAL
  [../]

[]

[AuxKernels]
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
  [./bnd_li_flux]
    type = DiffusionFluxNormalToBoundaryAux
    variable = bndliflux
    boundary = 'top'
    diffusion_variable = V
    diffusivity = thermal_conductivity
  [../]
  [./chemical_potential]
    type = ADMaterialRealAux
    variable = chemical_potential
    property = chemical_potential
  [../]

  [./equil_potential]
    type = ADMaterialRealAux
    variable = equil_potential
    property = equilibrium_potential
  [../]

[]

[Kernels]
  [./li_ion_V]
    type = ADHeatConduction
    variable = V
    use_displaced_mesh = false
  [../]
  [./Diffusion]
#    type = ADChemoMechanoAnsioDiffusion
    type = ADChemoMechanoDiffusion
    variable = li_metal_conc
#    diffusivity = diffusivity
    diffusivity = 3e-3
  [../]
  [./diffusion_dt]
    type = ADTimeDerivative
    variable = li_metal_conc
  [../]
[]

[Materials]
#  [./diffusivity_Li]
#    type = ADDiffusionAlongPrincipalDirectionsMaterial
#    diffusivity_vector = '3e-3 0 0'
#  [../]

  [./equilibrium_potential]
    type = ADComputeEquilibriumPotential
    R = 8.31446
    faraday = 96.4853329
    temperature = 298
    cref = 0.0528
    concentration = li_metal_conc
    include_conc = true
    include_reaction_rate = true
    # reaction_rate = 780.0
    reaction_rate_function = reaction_rate
    include_mechanical_effects = false
    exclude_elastic_contribution = true
  [../]

  [./thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 1
  [../]
[]
[Functions]
  # --- This is the OCV ---
  [./reaction_rate]
    type = PiecewiseLinear
    data_file = data/eqmV_AgLi.csv
    format = columns
    scale_factor = 1.0
    xy_in_file_only = false
    x_index_in_file = 0
    y_index_in_file = 2
  [../]
[]


[BCs]
  [./current]
    type = ADButlerVolmerBC
    current_density = 1.2e-3
    exchange_current_density = 1e-6
    faraday = 96.4853329
    Temperature = 298
    R = 8.314462681
    variable = V
    boundary = top
  [../]
  [./conc]
    type = ADNeumannBC
    variable = li_metal_conc
    value = 1.24371375e-5
    boundary = top
  [../]
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
    boundary = top
    variable = V
  [../]
  [./V0]
    type = SideAverageValue
    boundary = top
    variable = equil_potential
  [../]
  [./chem_potential]
    type = SideAverageValue
    boundary = top
    variable = chemical_potential
  [../]
  [./conc]
    type = SideAverageValue
    boundary = top
    variable = li_metal_conc
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  automatic_scaling = true
  compute_scaling_once = true
  petsc_options_iname = '-pc_type -pc_mat_solver_package -snes_linesearch_type -snes_force_iteration'# -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu superlu_dist basic 1'#     NONZERO               1e-20               '
  dt = 10
  nl_max_its = 35
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-6
  dtmax = 20
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 10
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 40
    # timestep_limiting_postprocessor = matl_ts_min
  [../]
  end_time = 8000

[]

[Outputs]
  exodus = true
  file_base = rst/test_ag
  csv = true
[]
