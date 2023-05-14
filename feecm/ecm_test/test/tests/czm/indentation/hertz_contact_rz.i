# Hertz Contact: Sphere on sphere

# Spheres have the same radius, Young's modulus, and Poisson's ratio.

# Define E:
# 1/E = (1-nu1^2)/E1 + (1-nu2^2)/E2
#
# Effective radius R:
# 1/R = 1/R1 + 1/R2
#
# F is the applied compressive load.
#
# Area of contact a::
# a^3 = 3FR/4E
#
# Depth of indentation d:
# d = a^2/R
#
#
# Let R1 = R2 = 2.  Then R = 1.
#
# Let nu1 = nu2 = 0.25, E1 = E2 = 1.40625e7.  Then E = 7.5e6.
#
# Let F = 10000.  Then a = 0.1, d = 0.01.
#

## Note: There is not a good way to check the result.  The standard approach is
## to map contact pressure as a function of radius, but we don't have the
## contact pressure available.  See the description on Wikipedia for details of
## analytic equations, and the Abaqus Benchmarks Manual, 1.1.11, for a plot of
## contact pressure vs. radius.
[GlobalParams]
  volumetric_locking_correction = true
  displacements = 'disp_x disp_y'
[]

[Problem]
  coord_type = RZ
[]

[Mesh] #Comment
  displacements = 'disp_x disp_y'
  allow_renumbering = false # Mesh
  [mesh]
    type = FileMeshGenerator
    file = 'hertz_contact_rz2.e'
  []
  [bottom]
    type = SideSetsFromNormalsGenerator
    input = 'mesh'
    normals = '0 -1 0'
    new_boundary = 'test'
  []
[]

[Functions]
  [pressure]
    type = PiecewiseLinear
    x = '0. 1. 2.'
    y = '0. 1. 1.'
    scale_factor = 795.77471545947674 # 10000/pi/2^2
  []
  [disp_y]
    type = PiecewiseLinear
    x = '0.  1.    2.'
    y = '0. -0.01 -0.01' # Functions
  []
[]

[Variables]

  [disp_x]
    order = FIRST
    family = LAGRANGE
  []

  [disp_y]
    order = FIRST
    family = LAGRANGE
  []
  # Variables
[]

[Modules/TensorMechanics/Master]
  [all]
    add_variables = true
    strain = FINITE
    block = '1 1000'
    use_automatic_differentiation = true
    generate_output = 'stress_xx stress_yy stress_zz vonmises_stress'
  []
[]

[BCs]

  [base_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'test'
    value = 0.0
  []

  [symm_x]
    type = DirichletBC
    variable = disp_x
    boundary = 1
    value = 0.0
  []
  [disp_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 2
    function = disp_y
  []
  # BCs
[]

[Contact]
  [mech]
    primary = 100
    secondary = 1000
    model = frictionless
    formulation = mortar
    mortar_approach = legacy
    # model = coulomb
    # formulation = penalty
    # normalize_penalty = true
    # friction_coefficient = 0.4
    # penalty = 8e7
    # tangential_tolerance = 0.005
  []
[]

[Materials]

  [tensor]
    type = ADComputeIsotropicElasticityTensor
    block = '1'
    youngs_modulus = 5000e3
    poissons_ratio = 0.33
  []
  [stress]
    type = ADComputeFiniteStrainElasticStress
    block = '1'
  []

  [tensor_1000]
    type = ADComputeIsotropicElasticityTensor
    block = '1000'
    youngs_modulus = 200e3
    poissons_ratio = 0.3
  []
  [stress_1000]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plas'
    block = '1000'
  []
  [plas]
    type = ADIsoTropicHyperVisco
    absolute_tolerance = 1e-6
    # block = 0
    # relative_tolerance = 1e-06
    hardening_exponent = 1.0
    saturation_resistance = 2.0e3
    initial_resistance = 1e3
    hardening_modulus = 100.0
    rate_exponent = 0.01
    reference_strain_rate = 0.1
    effective_inelastic_strain_name = effective_plastic_strain
    max_inelastic_increment = 1.0
    internal_solve_full_iteration_history = true
    internal_solve_output_on = on_error
    block = '1000'
  []
  # Materials
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]

  type = Transient

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'
  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type '
  petsc_options_value = 'lu'

  # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  # petsc_options_value = 'hypre    boomeramg      101'

  line_search = 'none'

  nl_abs_tol = 1e-10

  l_max_its = 200

  start_time = 0.0
  dt = 0.5
  end_time = 2.0 # Executioner
  automatic_scaling = true
  resid_vs_jac_scaling_param = 0.5
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-3
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 30
    # timestep_limiting_postprocessor = matl_ts_min
  []
[]

[Postprocessors]
  [maxdisp]
    type = NodalVariableValue
    nodeid = 39 # 40-1 where 40 is the exodus node number of the top-left node
    variable = disp_y
  []
[]

[Outputs]
  [out]
    type = Exodus
    elemental_as_nodal = true # Outputs
  []
[]
