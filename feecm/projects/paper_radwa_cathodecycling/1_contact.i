!include sync_times.i
## ============================================================================
## UNIT SYSTEM: um (length) - MPa (stress) - uN (force).  1 MPa = 1 uN/um^2.
##   Representative cell = 5 x 5 um; NMC particle radius ~3.9 um.  The mesh is
##   unit-agnostic (coords read 5.0); these material params set the um scale.
## Catholyte / solid electrolyte = LPSC (Li6PS5Cl).
##
## CYCLING: the NMC eigenstrain is driven by a triangular temperature wave of
##   period cycle_period, repeated n_cycles times (expansion on the first half
##   of each cycle, contraction on the second half).
## ============================================================================
## --- Catholyte (LPSC solid electrolyte) properties ---
## Rate-independent J2 plasticity: linear elastic + isotropic (von Mises) yield,
## yield strength tied to the measured Vickers hardness (Tabor sigma_y ~ H_v/3).
##
## Run: mpiexec -np 8 ./contact_loss-opt -i 1_contact.i
##
## Material Properties:
ymod_lpsc=550              # Young Modulus of LPSC [MPa] (22 GPa)
Hv_lpsc=2000               # Vickers Hardness of LPSC [MPa] (2 GPa)
pr_lpsc=0.37               # Poissons Ratio of LPSC

ymod_nmc=177500           # Young Modulus of NMC [MPa] (177.5 GPa, isotropic elastic cathode)
pr_nmc=0.33               # Poissons Ratio of NMC

## --- Calculated Material Properties ---
ustr_lpsc=${fparse Hv_lpsc * 0.012}                 # Ultimate Strength of SE [MPa] (Coe. = 83.33, what Shafee used)
ystr_lpsc=${fparse ustr_lpsc / 1.2}                 # Yield Strength of SE [MPa]
plstr=${fparse (ustr_lpsc - ystr_lpsc) / (1000)}    # Plastic Strain of SE (Plastic Hardening Tangent is 1000 MPa)

## --- Boundary Conditions Properties ---
sptop=7                       # Stack Pressure [MPa]
alpha_nmc=2.92245917e-4       # Thermal expansion coefficient of NMC (8.1% expansion)

## --- Cycling Control ---
n_cycles=1000                                   # Number of charge/discharge cycles
cycle_period=1.0                                # Duration of one full expand + contract cycle
temp_peak=90.0                                  # Temperature at full lithiation (end of expansion)
t_end=${fparse n_cycles * cycle_period}         # Total simulated time
dt_cycle=${fparse cycle_period / 40}            # 40 steps per cycle resolves the load reversal


[Problem]
  type = FEProblem
  solve = true
[]
[Mesh]
        [./fmg]
                type = FileMeshGenerator
                ## Mesh (AM:LPSC = 70/30 Weight ratio)
                file = mesh_7030.msh
        []
[]
[GlobalParams]
        displacements = 'disp_x disp_y'
[]
[AuxVariables]
        [./temp]
                initial_condition = 0
        [../]
        [./eigenstrain_xx]
                order = FIRST
                family = MONOMIAL
                block = 'block_NMC'
        [../]
        [./eigenstrain_yy]
                order = FIRST
                family = MONOMIAL
                block = 'block_NMC'
        [../]
        [./total_strain_xx]
                order = FIRST
                family = MONOMIAL
                block = 'block_NMC'
        [../]
        [./total_strain_yy]
                order = FIRST
                family = MONOMIAL
                block = 'block_NMC'
        [../]
[]
[Functions]
        [./temperature_load]
                type = ParsedFunction
                # Triangular wave: 0 -> peak at half period -> 0, repeated every period.
                expression = '2*peak/per * (per/2 - abs(t - per*floor(t/per) - per/2))'
                symbol_names = 'peak per'
                symbol_values = '${temp_peak} ${cycle_period}'
        [../]
        [./press_ramp]
                type = PiecewiseLinear
                x = '0    0.05  ${t_end}'
                y = '0.0  1.0   1.0'
        [../]
        [./hf_LPSC]
                type = PiecewiseLinear
                x = '0 ${plstr}'
                # Hardening Function LPSC
                y = '${ystr_lpsc} ${ustr_lpsc}'
        [../]
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [./NMC]
        strain = FINITE
        decomposition_method = EigenSolution
        add_variables = true
        eigenstrain_names = eigenstrain
        generate_output = 'stress_xx stress_yy vonmises_stress strain_xx strain_yy'
        block = 'block_NMC'
      [../]
      [./LPSC]
        strain = FINITE
        decomposition_method = EigenSolution
        add_variables = true
        generate_output = 'stress_yy stress_xx vonmises_stress plastic_strain_xx plastic_strain_yy'
        block = 'block_LPSC'
      [../]
    [../]
  [../]
[]

[AuxKernels]
        [./tempfuncaux]
                type = FunctionAux
                variable = temp
                function = temperature_load
        [../]
        [./eigenstrain_yy]
                type = RankTwoAux
                block = 'block_NMC'
                rank_two_tensor = eigenstrain
                variable = eigenstrain_yy
                index_i = 1
                index_j = 1
                execute_on = 'initial timestep_end'
        [../]
        [./eigenstrain_xx]
                type = RankTwoAux
                block = 'block_NMC'
                rank_two_tensor = eigenstrain
                variable = eigenstrain_xx
                index_i = 0
                index_j = 0
                execute_on = 'initial timestep_end'
        [../]
        [./total_strain_yy]
                type = RankTwoAux
                block = 'block_NMC'
                rank_two_tensor = total_strain
                variable = total_strain_yy
                index_i = 1
                index_j = 1
                execute_on = 'initial timestep_end'
        [../]
        [./total_strain_xx]
                type = RankTwoAux
                block = 'block_NMC'
                rank_two_tensor = total_strain
                variable = total_strain_xx
                index_i = 0
                index_j = 0
                execute_on = 'initial timestep_end'
        [../]
[]

[Contact]
        [nmc_lpsc]
                primary = 'block_NMC_right'
                secondary = 'block_LPSC_left'
                penalty = ${ymod_lpsc}
                formulation = penalty
                tangential_tolerance = 0.0001
        []
[]

[BCs]
        [./x_disp]
                type = DirichletBC
                variable = disp_x
                boundary = 'block_left block_right'
                value = 0.0
        [../]
        [./y_disp]
                type = DirichletBC
                variable = disp_y
                boundary = 'block_bottom'
                value = 0.0
        [../]
        [./top_press]
                type = Pressure
                variable = disp_y
                boundary = 'block_top'
                factor   = ${sptop}
                function = press_ramp
        [../]
[]

[Materials]
        [./elasticity_tensor_NMC]
                type = ComputeIsotropicElasticityTensor
                block = 'block_NMC'
                youngs_modulus = ${ymod_nmc}
                poissons_ratio = ${pr_nmc}
        [../]
        [./elastic_stress_NMC]
                # NMC cathode is modeled as isotropic linear elastic (no plasticity).
                type = ComputeFiniteStrainElasticStress
                block = 'block_NMC'
        [../]
        [./thermal_expansion_strain_NMC]
                type = ComputeThermalExpansionEigenstrain
                block = 'block_NMC'
                stress_free_temperature = 0
                thermal_expansion_coeff = ${alpha_nmc}
                temperature = temp
                eigenstrain_name = eigenstrain
        [../]
        [./elasticity_tensor_LPSC]
                type = ComputeIsotropicElasticityTensor
                block = 'block_LPSC'
                youngs_modulus = ${ymod_lpsc}
                poissons_ratio = ${pr_lpsc}
        [../]
        [./isotropic_plasticity_LPSC]
                type = IsotropicPlasticityStressUpdate
                block = 'block_LPSC'
                yield_stress = ${ystr_lpsc}
                hardening_function = hf_LPSC
        [../]
        [./radial_return_stress_LPSC]
                type = ComputeMultipleInelasticStress
                tangent_operator = nonlinear
                inelastic_models = 'isotropic_plasticity_LPSC'
                block = 'block_LPSC'
        [../]
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  automatic_scaling = true
  solve_type = NEWTON
  dtmin = 1e-6
  dtmax = ${dt_cycle}
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type -mat_mumps_icntl_24 -mat_mumps_icntl_14'
  petsc_options_value = 'lu       mumps                       1                           200'
  line_search = none
  nl_max_its = 99
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-7
  l_tol = 1e-8
  start_time = 0.0
  n_startup_steps = 1
  end_time = ${t_end}
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    optimal_iterations = 8
    iteration_window  = 2
    growth_factor = 2.0
    cutback_factor = 0.7
    cutback_factor_at_failure = 0.5
    linear_iteration_ratio = 100
  []
[]

[Outputs]
  file_base = rst/E${ymod_lpsc}_H${Hv_lpsc}_spTop${sptop}_N${n_cycles}
  checkpoint = true
  [exodus]
    type = Exodus
    sync_times = ${output_times_str}
    sync_only = true
  []
[]