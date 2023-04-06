# Initial test for elastic cathode material
# No mechanics are included in this
# Reaction rate obtained from comsol NMC electrode
[Mesh]
  [./Si]
    type = GeneratedMeshGenerator
    nx = 100
    ny = 2
    xmin = 0
    xmax = 500
    ymin = 0
    ymax = 0.01
    elem_type = QUAD4
    dim = 2
  [../]
[]
# [GlobalParams]
#   displacements = 'ux uy'
# []

[Variables]
  [./V]
  [../]
  [./li_metal_conc]
    initial_condition = 614.172e-6
  [../]

[]

[Kernels]
  [./li_ion_V]
    type = ADHeatConduction
    variable = V
    use_displaced_mesh = false
  [../]
  [./Diffusion]
    type = ADMatDiffusion
    variable = li_metal_conc
    # stress_based_chemical_potential = mu_sigma
    diffusivity = diffusivity
    # R = 96.4853329
    # temperature = 298
  [../]
  [./diffusion_dt]
    type = ADTimeDerivative
    variable = li_metal_conc
  [../]
[]

[Materials]

  [./diffusivity_Li]
    type = ADGenericConstantMaterial
    prop_names = 'diffusivity'
    prop_values = '3e-2'
  [../]

  [./equilibrium_potential]
    type = ADComputeEquilibriumPotential
    R = 8.31446
    faraday = 96.4853329
    temperature = 298
    cref = 7.874e-2
    concentration = li_metal_conc
    include_conc = false
    include_reaction_rate = true
    # reaction_rate = 780.0
    reaction_rate_function = reaction_rate
    include_mechanical_effects = false
    exclude_elastic_contribution = true
  [../]


  [./thermal_conductivity1]
    type = ADHeatConductionMaterial
    thermal_conductivity = 10
  [../]
[]

[Functions]
  [./reaction_rate]
    type = PiecewiseLinear
    format = columns
    data_file = 'nmc_equilibrium_potential.csv'
  [../]
[]


[BCs]
  [./current]
    type = ADButlerVolmerBC
    current_density = 1
    exchange_current_density = 1e-6
    faraday = 96.4853329
    Temperature = 298
    R = 8.314462681
    # equilibrium_potential = V0
    variable = V
    boundary = top
    # extra_vector_tags = 'ref'
  [../]
  [./conc]
    type = ADNeumannBC
    variable = li_metal_conc
    # v = bndliflux
    value = 1.24371375e-5
    boundary = top
    # extra_vector_tags = 'ref'
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
  dt = 200
  # num_steps = 20
  nl_max_its = 35
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-6
  dtmax = 200
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
