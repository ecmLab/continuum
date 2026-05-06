# Simplified Phase 2 EDL benchmark without AD type mismatches

[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0.0
  xmax = 3.04e-8          # 10 Debye lengths (~30 nm)
  nx = 100
  # Remove bias for now to avoid degenerate elements
[]

[Variables]
  [./phi]
    initial_condition = 0.0
    scaling = 1.0e-06    # Scale down phi significantly to balance with charge source
  [../]
  [./c_plus]
    initial_condition = 10.0   # mol/m^3
    scaling = 1.0e+05    # Scale up concentrations to balance diffusion terms
  [../]
  [./c_minus]
    initial_condition = 10.0
    scaling = 1.0e+05    # Scale up concentrations to balance diffusion terms
  [../]
[]

[Materials]
  [./constants]
    type = GenericConstantMaterial
    prop_names = 'D_plus D_minus permittivity F R T z_plus z_minus'
    prop_values = '1.0e-9 1.0e-9 6.95104e-10 96485.33212 8.314462618 298.15 1.0 -1.0'
  [../]
[]

[Kernels]
  # Poisson equation: -∇·(ε∇φ) = F (z_+ c_+ + z_- c_-)
  [./phi_diffusion]
    type = MatDiffusion
    variable = phi
    diffusivity = permittivity
  [../]
  [./phi_source_plus]
    type = CoupledForce
    variable = phi
    v = c_plus
    coef = -96485.33212
  [../]
  [./phi_source_minus]
    type = CoupledForce
    variable = phi
    v = c_minus
    coef = 96485.33212
  [../]

  # Steady Nernst–Planck for c_+
  [./cplus_diffusion]
    type = MatDiffusion
    variable = c_plus
    diffusivity = D_plus
  [../]
  [./cplus_electromigration]
    type = CoupledForce
    variable = c_plus
    v = phi
    coef = -9.649e-13  # Approximated electromigration coefficient
  [../]

  # Steady Nernst–Planck for c_-
  [./cminus_diffusion]
    type = MatDiffusion
    variable = c_minus
    diffusivity = D_minus
  [../]
  [./cminus_electromigration]
    type = CoupledForce
    variable = c_minus
    v = phi
    coef = 9.649e-13   # Approximated electromigration coefficient (opposite sign)
  [../]
[]

[BCs]
  # Simplified Robin BC at electrode (left)
  [./phi_left]
    type = DirichletBC
    variable = phi
    boundary = left
    value = 1.0e-3  # 1 mV applied potential
  [../]

  # Reference gauge in bulk
  [./phi_right]
    type = DirichletBC
    variable = phi
    boundary = right
    value = 0.0
  [../]

  # Bulk concentrations at right boundary
  [./cplus_right]
    type = DirichletBC
    variable = c_plus
    boundary = right
    value = 10.0
  [../]
  [./cminus_right]
    type = DirichletBC
    variable = c_minus
    boundary = right
    value = 10.0
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  nl_rel_tol = 1e-5    # Relax tolerance for demonstration
  nl_abs_tol = 5e-9    # Set absolute tolerance slightly above balanced residual scale
  nl_max_its = 100     # Increase iteration limit
  l_tol = 1e-10
  l_max_its = 400
  petsc_options_iname = '-pc_type -ksp_type'
  petsc_options_value = 'lu       preonly'
[]

[Debug]
  show_var_residual_norms = true
[]

[Outputs]
  exodus = true
  csv = true
  file_base = rst/comsol_phase2_edl_simplified
[]