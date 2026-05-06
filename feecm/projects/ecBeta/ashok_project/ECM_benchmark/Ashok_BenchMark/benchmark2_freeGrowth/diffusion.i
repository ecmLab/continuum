[Mesh]
    patch_size = 80
    patch_update_strategy = auto
    parallel_type = REPLICATED
    [./mesh]
      type = GeneratedMeshGenerator
      xmax = 5.0
      xmin = 0.0
      ymax = 1.0
      ymin = 0.0
      dim = 2
      elem_type = QUAD4
      nx = 5
      ny = 2
    [../]
  []
  
  
  
  [Problem]
    coord_type = 'RZ'
  []
  [Variables]
  
  
    [./conc]
     # initial_condition = 0
  #    block = 'interLayer'
    [../]
  []
  
  [AuxVariables]
    [flux_x]
      order = CONSTANT
      family = MONOMIAL
  #    block = interLayer
    []
    [flux_y]
      order = CONSTANT
      family = MONOMIAL
  #    block = interLayer
    []
  
  []
  
 
  
  [AuxKernels]
    [./flux_x]
      type = AnisoTropicDiffusionFluxAux
      variable = flux_x
      diffusion_variable = conc
      component = x
    [../]
    [./flux_y]
      type = AnisoTropicDiffusionFluxAux
      variable = flux_y
      diffusion_variable = conc
      component = y
    [../]
  []
  
  [Kernels]
    [./diffusion]
      type = ADMatAnisoDiffusion
      variable = conc
      diffusivity = diffusivity
      use_displaced_mesh = false
  #    block = 'interLayer'
    [../]
    [./li_metal_dt]
      type = ADTimeDerivative
      variable = conc
      use_displaced_mesh = false
  #    block = 'interLayer'
    [../]
  
  []
  
  [Materials]
    [./diffusivity_Li1]
      type = ADDiffusionAlongPrincipalDirectionsMaterial
  #   block = 'interLayer'
      diffusivity_vector = '0 1e8 0'
    [../]
    [swell_normal]
        type = ADGenericConstantVectorMaterial
        #boundary = 'bottom'
        prop_names = swell_normal
        prop_values = '0 1 0'
    []
  []
  
  [BCs]
    # [./top_flux]
    #   type = ADNeumannBC
    #   variable = conc
    #   value = 1
    #   boundary = bottom
    # [../]
    # [./top_bc]
    #   type = ADDirichletBC
    #   variable = conc
    #   value = 0
    #   boundary = top
    # [../]
    [./top_bc]
      type = ADDirichletBC
      variable = conc
      value = 0
      boundary = top
    [../]
      [./bottom_bc]
        type = ADNeumannBC
        variable = conc
        value = 1
        boundary = bottom
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
    petsc_options = '-snes_converged_reason -ksp_converged_reason'
    petsc_options_iname = '-pc_type -pc_mat_solver_package'
    petsc_options_value = 'lu superlu_dist'
    l_max_its = 50
    nl_max_its = 25
  
    nl_rel_tol = 1e-10
    nl_abs_tol = 1e-8
  
    start_time = 0.0
    dt = 1
    dtmax = 2.0
    dtmin = 1e-9
    #num_steps = 1000
    end_time = 1000
  
  [] # Executioner
  
  [Outputs]
    [./out]
      type = Exodus
      file_base = rst/diffusion_dummy
    [../]
  [] # Outputs
  
  [Debug]
    show_var_residual_norms = true
    show_material_props = true
  []
  