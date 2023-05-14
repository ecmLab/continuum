# Simulation of Li and electron transport under the influence of voltage difference
# Variables are Li and electron concentration (cLi, cEl) and potential distribution (phi)
# Diffusivity of Li-ions in LLZO is D_Li = 5*10^5 nm^2/s, Diffusivity of electron in LLZO is D_e = 50 nm^2/s
# Permissivity of LLZO is epsilon_b = 50*epsilon_0 = 50*8.85*10^-12 F/m = 4.43*10^-19 F/nm
# A constant epsilong_e is defined with permissivity divided by the electron charge: epsilon_e = epsilong_b/e = 2.7613, in unit 1/(V*nm)
# Length in unit nm

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
  xmax = 100 #Length of the LLZO material, in unit nm
  ymax = 100 #height of the LLZO material, in unit nm
[]

[Variables]
  [./cE]
  initial_condition = 0.0 # in unit (#/nm^3)
#  initial_condition = 1  # bulk concentration of Li in LLZO ~0.167*10^4 mol/m^3 = 6 (#/nm^3)
  [../]
  [./phi]
  initial_condition = 0.0
  [../]
[]

[Functions]
  active = 'source_func'
  [./source_func]
    type = PiecewiseMultilinear
    data_file = chrgDen.txt
#    type = ParsedFunction
#    value = 'if(x<5.5,1,0)'
  [../]
[]

[Kernels]
# This kernel describes the Laplace term in the Nernst-Planck equation
# diffusivity is electron diffusivity in llzo
  [./electron_diffusion]
#    type = ParamDiffusion
    type = MatDiffusion
    variable = cE
    diffusivity = diff_coefficient
  [../]
  [./dcE_dt]
    type = TimeDerivative
    variable = cE
  [../]

# This kernel describes the convection term in the Nernst-Planck equation
# diffusivity is electron diffusivity in llzo
  [./electrical_convection]
    type = NernstPlanckConvection
    variable = cE
    diffusivity = diff_coefficient
    Voltage = phi
    zIons    = -1
  [../]

# This kernel describes the convection term in the Nernst-Planck equation
# diffusivity is ratio of electron diffusivity and the normalized permittivity
  [./electron_source]
    type = sourceTerm
    variable = cE
    srcCoef  = src_coefficient
    function = source_func
    block = 0
   # function = volumetric_c
  [../]

# This kernel describes the Laplace term in the Poisson equation
# diffusivity is the permittivity in LLZO normolized by electron charge: epsilon_e = epsilong_b/e, in unit 1/(V*nm)
# scale is the scale up of this term, for problem stability
  [./Voltage_Poisson]
    type = MatDiffusion
    variable = phi
    diffusivity = permittivity_e
    scale    = 1
  [../]

# This kernel describes the charge density in the Poisson equation caused by the transport of electron concentration
# zIons is the charge state, zIons = -1 for electron
# scale is the scale up of this term, for problem stability
  [./Voltage_chargeDensity]
    type = chargeDensity
    variable = phi
    conIons  = cE
    zIons    = -1
    scale    = 1
  [../]

[]

[BCs]
#  [./concentration_bottom]
#    type = FunctionDirichletBC
    #type = NeumannBC
#    variable = cE
#    boundary = right
#    #value = 1.5e-9 #1.76e-10 #2.77e-6 #wt %/s (mol/cc^3/s) kd=8.2e-6 m/s del c =6.6e-3 mol/cc sgb =1.14e-2
#    #value = 0.0  #wt %
#    function = bc_func
#  [../]
#  [./concentration_bottom]
#    type = DirichletBC
#    variable = cE
#    boundary = 'bottom'
#    value    = 0.0
#  [../]
#  [./flux_top]
#    type = NeumannBC
#    variable = cE
#    boundary = top
#    value = 100 # (#/nm^3)
#  [../]
  [./voltage_bottom]
   type = DirichletBC
    variable = phi
    boundary = bottom
    value = 0.0 #520.0 # (K) 
  [../]
#  [./voltage_left]
#   type = DirichletBC
#    variable = phi
#    boundary = left
#    value = 0.7e-04 #520.0 # (K) 
#  [../]
[]

[VectorPostprocessors]
  [left_u]
    type = SideValueSampler
    variable = cE
    boundary = left
    sort_by = y
  []
[]

[Materials]
  [./diff_coeff]
## Unit of parameters:
# Diffusivity of Li-ions in LLZO is D_Li = 5*10^5 nm^2/s, Diffusivity of electron in LLZO is D_e = 50 nm^2/s
# permissivity_e is defined with permissivity divided by the electron charge: epsilon_e = epsilong_b/e, in unit 1/(V*nm)
#in m^2/s mol/cc 3.2e-9 6.8e-3 6.44e-9 mobility = mob_coefficient= ( Dze/(k_b*T))
   type = GenericFunctionMaterial
   prop_names = 'diff_coefficient  conc_sat  src_coefficient  permittivity_e'
   prop_values = '50               6.8e-3    18.11            2.7613' 
   block = 0
  [../]
 []

[Preconditioning]
 [./coupled]
    type = SMP
    full = true
  [../]
[]


[Executioner]
  type = Transient
  num_steps = 500
  dt = 0.001
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart'
  petsc_options_value = 'asm lu 1 101'
#  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
#  petsc_options_value = 'hypre    boomeramg      101'
[]

[Outputs]
  exodus = true
  csv = true
  file_base = rst/out
[]
