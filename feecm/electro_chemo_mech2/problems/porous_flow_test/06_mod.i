# Darcy flow with a tracer single material
[Mesh]
  [block]
    type = GeneratedMeshGenerator
    xmin = 0
    xmax = 10
    dim = 2
    ymax = 10
    ymin = 0
    elem_type = QUAD4
    bias_x = 1.0
    nx = 10
    ny = 10
    # nr = 10
    # rmin = 1.0
    # rmax = 10
    # growth_r = 1.4
    # nt = 4
    # dmin = 0
    # dmax = 90
  []
  [aquifer]
    type = SubdomainBoundingBoxGenerator
    block_id = 1
    bottom_left = '0 -0.1 0'
    top_right = '5 20 0'
    input = block
  []
  [rename]
    type = RenameBlockGenerator
    old_block_id = '0 1'
    new_block_name = 'caps aquifer'
    input = 'aquifer'
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [porepressure]
  []
  [tracer_concentration]
  []
[]

[ICs]
  [tracer_concentration]
    type = FunctionIC
    # function = '0.5*if(x < 1e-6, if (y <= 5.0, 1.0, 0.0),0)'
    function = '0.5*if(x < 1e-6, 1,0)'
    # function = '0.5*if (x < 1e-6, if (y < 5,1, 0))'
    variable = tracer_concentration
  []
[]

[PorousFlowFullySaturated]
  porepressure = porepressure
  coupling_type = Hydro
  gravity = '0 0 0'
  fp = the_simple_fluid
  mass_fraction_vars = tracer_concentration
  stabilization = KT # Note to reader: 06_KT.i uses KT stabilization - compare the results
  flux_limiter_time = superbee
[]

[BCs]
  # [constant_injection_porepressure]
  #   type = DirichletBC
  #   variable = porepressure
  #   value = 1E6
  #   boundary = left
  #   # boundary = bottom
  #   # function = '1e6*if (x <= 5, 1, 0)'
  # []
  [constant_outer_porepressure]
    type = DirichletBC
    variable = porepressure
    value = 0
    boundary = right
  []
  [constant_flux_injection]
    type = PorousFlowSink
    boundary = left
    variable = tracer_concentration
    flux_function = 0.1
  []
  # [injected_tracer]
  #   type = FunctionDirichletBC
  #   variable = tracer_concentration
  #   boundary = left
  #   # value =
  #   # function = '0.5*if(x < 1e-6, if (y <= 5.0, 1.0, 0.0),0)'
  #   function = '0.5*if (x <= 1e-6, 1, 0)'
  # []
[]

[Modules]
  [FluidProperties]
    [the_simple_fluid]
      type = SimpleFluidProperties
      bulk_modulus = 2E9
      viscosity = 1.0E-3
      density0 = 1000.0
      
    []
  []
[]

[Materials]
  [porosity_aquifer]
    type = PorousFlowPorosityConst
    porosity = 0.3
    block = 'aquifer'
    # mechanical = false
    # fluid = true
  []
  [porosity_caps]
    type = PorousFlowPorosityConst
    porosity= 0.1
    block = 'caps'
    # mechanical = false
    # fluid = false
  []

  [permeability_aquifer]
    type = PorousFlowPermeabilityConst
    block = aquifer
    permeability = '1E-12 0 0   0 1E-12 0   0 0 1E-12'
  []
  [permeability_caps]
    type = PorousFlowPermeabilityConst
    block = caps
    permeability = '1E-15 0 0   0 1E-15 0   0 0 1E-15'
  []
[]

[Preconditioning]
  active = preferred_but_might_not_be_installed
  [basic]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = ' asm      lu           NONZERO                   2'
  []
  [preferred_but_might_not_be_installed]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       superlu_dist'
  []
[]

[Executioner]

  type = Transient
  automatic_scaling = true
  solve_type = Newton
  end_time = 3E6
  dt = 1E5
  nl_rel_tol = 1E-11
[]

[Outputs]
  exodus = true
[]
