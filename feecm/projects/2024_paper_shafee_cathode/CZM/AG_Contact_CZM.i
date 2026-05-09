## ----------------------------------------------------------------------------
##  AG_Contact_CZM.i
##  NMC (swelling) / LPS contact problem with a Bilinear Mixed-Mode CZM
##  replacing the original penalty contact formulation.
## ----------------------------------------------------------------------------

## Bulk material parameters (unchanged from AG_Contact_Loss.i) ---------------
## Young modulus of LPSCL [MPa]
ymod = 900
## Yield strength of LPSCL [MPa]
ystr = 196
ystr_plus_4 = ${fparse ystr + 4.0}
## Stack pressure [MPa]
sptop = 6

## Cohesive zone parameters --------------------------------------------------
## Peak normal traction (cohesive strength, Mode I)         [MPa]
sigma_max = 20
## Peak shear traction  (cohesive strength, Mode II)        [MPa]
tau_max   = 20
## Mode-I  fracture energy release rate                     [N/mm]
GIc       = 0.1
## Mode-II fracture energy release rate                     [N/mm]
GIIc      = 0.1
## Initial elastic (penalty) stiffness of the bonded zone   [MPa/mm]
## Rule of thumb: K * h_elem >> E_bulk so that pre-damage compliance is
## negligible.  1e5 MPa/mm gives a stiff bond for h_elem ~ 1e-3 mm.
Kpen      = 1e5
## Benzeggagh-Kenane mixed-mode exponent (1.0 - 2.0 typical)
eta_BK    = 2.0
## Numerical viscosity for damage regularisation (small, helps convergence)
visc_czm  = 1e-5

[Problem]
  type  = FEProblem
  solve = true
[]

[Mesh]
  ## The Gmsh file (input_mesh_Shafee.msh) defines block_NMC and block_LPS
  ## as two distinct physical surfaces.  Their interface is represented by
  ## TWO different curves -- block_NMC_right (curve tag 1) on the NMC side
  ## and block_LPS_left (curve tag 9) on the LPS side -- which are
  ## geometrically coincident but topologically separate (different node
  ## ids).  The CZM InterfaceKernel needs a sideset whose two sides are
  ## reachable through shared mesh topology, so we must first fuse the
  ## duplicate interface nodes.
  ##
  ## Workflow:
  ##   1. fmg     : load the original mesh as-is.
  ##   2. mesh_NMC: keep only block_NMC by deleting block_LPS.
  ##   3. mesh_LPS: keep only block_LPS by deleting block_NMC.
  ##   4. stitched: stitch the two meshes, pairing the coincident
  ##                boundaries 'block_NMC_right' <-> 'block_LPS_left'.
  ##                StitchedMeshGenerator merges the duplicate nodes by
  ##                spatial matching, so the resulting mesh is fully
  ##                connected across the (now interior) interface.
  ##   5. break_interface: split the interior face between the two blocks
  ##                back into two element-faces with split nodes, and
  ##                expose them as a single sideset called 'interface'
  ##                (which the CZM material/action below references).
  [./fmg]
    type = FileMeshGenerator
    file = input_mesh_Shafee.msh
  [../]
  [./mesh_NMC]
    type  = BlockDeletionGenerator
    input = fmg
    block = 'block_LPS'
  [../]
  [./mesh_LPS]
    type  = BlockDeletionGenerator
    input = fmg
    block = 'block_NMC'
  [../]
  [./stitched]
    type                    = StitchedMeshGenerator
    inputs                  = 'mesh_NMC mesh_LPS'
    stitch_boundaries_pairs = 'block_NMC_right block_LPS_left'
    ## Outer boundaries (block_left, block_right, block_top, block_bottom)
    ## that exist on BOTH input meshes are concatenated under their shared
    ## name -- no extra plumbing needed for the BCs.
  [../]
  [./break_interface]
    type            = BreakMeshByBlockGenerator
    input           = stitched
    split_interface = false
    interface_name  = 'interface'
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
    order  = FIRST
    family = MONOMIAL
    block  = 'block_NMC'
  [../]
  [./eigenstrain_yy]
    order  = FIRST
    family = MONOMIAL
    block  = 'block_NMC'
  [../]
  [./total_strain_xx]
    order  = FIRST
    family = MONOMIAL
    block  = 'block_NMC'
  [../]
  [./total_strain_yy]
    order  = FIRST
    family = MONOMIAL
    block  = 'block_NMC'
  [../]
  ## Aux variables to visualise the CZM state on the interface ---------------
  [./czm_damage]
    order  = CONSTANT
    family = MONOMIAL
  [../]
  [./czm_traction_n]
    order  = CONSTANT
    family = MONOMIAL
  [../]
  [./czm_jump_n]
    order  = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]
  [./temperature_load]
    type       = ParsedFunction
    ## Loading-unloading ramp (same as original input)
    expression = 'if(t<=0.5, 180*t, (-1)*(t-0.5)*180 + 90.0)'
  [../]
  [./hf_NMC]
    type = PiecewiseLinear
    x    = '0'
    y    = '12600'
  [../]
  [./hf_LPS]
    type = PiecewiseLinear
    x    = '0 0.004'
    y    = '${ystr} ${ystr_plus_4}'
  [../]
[]

[Physics]
  [./SolidMechanics]
    ## Bulk volumetric mechanics ---------------------------------------------
    [./QuasiStatic]
      [./NMC]
        strain             = FINITE
        add_variables      = true
        eigenstrain_names  = eigenstrain
        generate_output    = 'stress_xx stress_yy vonmises_stress strain_xx strain_yy'
        block              = 'block_NMC'
      [../]
      [./LPS]
        strain          = FINITE
        add_variables   = true
        generate_output = 'stress_yy stress_xx vonmises_stress plastic_strain_xx plastic_strain_yy'
        block           = 'block_LPS'
      [../]
    [../]
    ## Cohesive zone (interface mechanics) -----------------------------------
    ## This action adds the CZMInterfaceKernel(s) for each displacement
    ## component on the 'interface' sideset.  The traction is computed by
    ## the BiLinearMixedModeTraction material below.
    [./CohesiveZone]
      [./czm_ik]
        boundary = 'interface'
        strain   = FINITE
      [../]
    [../]
  [../]
[]

[AuxKernels]
  [./tempfuncaux]
    type     = FunctionAux
    variable = temp
    function = temperature_load
  [../]
  [./eigenstrain_yy]
    type            = RankTwoAux
    block           = 'block_NMC'
    rank_two_tensor = eigenstrain
    variable        = eigenstrain_yy
    index_i         = 1
    index_j         = 1
    execute_on      = 'initial timestep_end'
  [../]
  [./eigenstrain_xx]
    type            = RankTwoAux
    block           = 'block_NMC'
    rank_two_tensor = eigenstrain
    variable        = eigenstrain_xx
    index_i         = 0
    index_j         = 0
    execute_on      = 'initial timestep_end'
  [../]
  [./total_strain_yy]
    type            = RankTwoAux
    block           = 'block_NMC'
    rank_two_tensor = total_strain
    variable        = total_strain_yy
    index_i         = 1
    index_j         = 1
    execute_on      = 'initial timestep_end'
  [../]
  [./total_strain_xx]
    type            = RankTwoAux
    block           = 'block_NMC'
    rank_two_tensor = total_strain
    variable        = total_strain_xx
    index_i         = 0
    index_j         = 0
    execute_on      = 'initial timestep_end'
  [../]
  ## Pull damage / traction / jump from the CZM material -------------------
  [./damage_aux]
    type           = MaterialRealAux
    variable       = czm_damage
    property       = damage
    boundary       = 'interface'
    execute_on     = 'initial timestep_end'
  [../]
  [./traction_n_aux]
    type           = MaterialRealVectorValueAux
    variable       = czm_traction_n
    property       = traction_global
    component      = 0
    boundary       = 'interface'
    execute_on     = 'initial timestep_end'
  [../]
  [./jump_n_aux]
    type           = MaterialRealVectorValueAux
    variable       = czm_jump_n
    property       = displacement_jump_global
    component      = 0
    boundary       = 'interface'
    execute_on     = 'initial timestep_end'
  [../]
[]

[BCs]
  [./x_disp]
    type     = DirichletBC
    variable = disp_x
    boundary = 'block_left block_right'
    value    = 0.0
  [../]
  [./y_disp]
    type     = DirichletBC
    variable = disp_y
    boundary = 'block_bottom'
    value    = 0.0
  [../]
  [./top_press]
    type     = Pressure
    variable = disp_y
    boundary = 'block_top'
    factor   = ${sptop}
  [../]
[]

[Materials]
  ## --- NMC bulk -----------------------------------------------------------
  [./elasticity_tensor_NMC]
    type           = ComputeIsotropicElasticityTensor
    block          = 'block_NMC'
    youngs_modulus = 177500
    poissons_ratio = 0.33
  [../]
  [./isotropic_plasticity_NMC]
    type               = IsotropicPlasticityStressUpdate
    block              = 'block_NMC'
    yield_stress       = 12600.0
    hardening_function = hf_NMC
  [../]
  [./radial_return_stress_NMC]
    type             = ComputeMultipleInelasticStress
    tangent_operator = elastic
    inelastic_models = 'isotropic_plasticity_NMC'
    block            = 'block_NMC'
  [../]
  [./thermal_expansion_strain_NMC]
    type                    = ComputeThermalExpansionEigenstrain
    block                   = 'block_NMC'
    stress_free_temperature = 0
    thermal_expansion_coeff = 0.00028
    temperature             = temp
    eigenstrain_name        = eigenstrain
  [../]

  ## --- LPS bulk -----------------------------------------------------------
  [./elasticity_tensor_LPS]
    type           = ComputeIsotropicElasticityTensor
    block          = 'block_LPS'
    youngs_modulus = ${ymod}
    poissons_ratio = 0.23
  [../]
  [./isotropic_plasticity_LPS]
    type               = IsotropicPlasticityStressUpdate
    block              = 'block_LPS'
    yield_stress       = ${ystr}
    hardening_function = hf_LPS
  [../]
  [./radial_return_stress_LPS]
    type             = ComputeMultipleInelasticStress
    tangent_operator = elastic
    inelastic_models = 'isotropic_plasticity_LPS'
    block            = 'block_LPS'
  [../]

  ## --- Cohesive zone traction-separation law ------------------------------
  ## Bilinear mixed-mode (BK) traction-separation:
  ##   * Linear elastic up to a quadratic stress envelope.
  ##   * Linear softening down to zero traction; total dissipated energy
  ##     equals the BK mixed-mode fracture energy.
  ##   * Compression handled by full penalty -> no separate [Contact] needed.
  [./czm_law]
    type              = BiLinearMixedModeTraction
    boundary          = 'interface'
    normal_strength   = ${sigma_max}
    shear_strength    = ${tau_max}
    GI_c              = ${GIc}
    GII_c             = ${GIIc}
    penalty_stiffness = ${Kpen}
    eta               = ${eta_BK}
    viscosity         = ${visc_czm}
    displacements     = 'disp_x disp_y'
  [../]
[]

[Executioner]
  type                 = Transient
  automatic_scaling    = true
  solve_type           = NEWTON
  petsc_options_iname  = '-pc_type'
  petsc_options_value  = 'lu'
  ## A bit of damping helps once the CZM enters the softening branch.
  line_search          = bt
  nl_max_its           = 99
  nl_rel_tol           = 1e-6
  nl_abs_tol           = 1e-9
  l_tol                = 1e-8
  start_time           = 0.0
  n_startup_steps      = 1
  end_time             = 1
  [TimeStepper]
    type                      = IterationAdaptiveDT
    dt                        = 1e-2
    optimal_iterations        = 10
    growth_factor             = 1.5
    cutback_factor            = 0.5
    cutback_factor_at_failure = 0.2
    linear_iteration_ratio    = 100
  []
[]

[Outputs]
  exodus    = true
  file_base = rst/czm_ystr${ystr}_E${ymod}_spTop${sptop}_smax${sigma_max}_GIc${GIc}
  [./csv]
    type       = CSV
    execute_on = 'final'
  [../]
[]

[Postprocessors]
  [./Gap]
    type     = PointValue
    point    = '0.0 0.00250001 0.0'
    variable = disp_y
  [../]
  [./Temp]
    type     = ElementAverageValue
    variable = temp
  [../]
  ## Maximum interface damage (1.0 = fully debonded somewhere on the interface)
  [./MaxDamage]
    type     = ElementExtremeValue
    variable = czm_damage
    boundary = 'interface'
  [../]
  ## Average normal traction across the interface
  [./AvgNormalTraction]
    type     = SideAverageValue
    variable = czm_traction_n
    boundary = 'interface'
  [../]
[]

[VectorPostprocessors]
  [./disp_xy_along_arc]
    type       = NodalValueSampler
    variable   = 'disp_x disp_y'
    boundary   = 'interface'
    sort_by    = x
    execute_on = 'final'
  [../]
[]