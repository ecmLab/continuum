## ============================================================================
## UNIT SYSTEM: um (length) - MPa (stress) - uN (force).  1 MPa = 1 uN/um^2.
##   Representative cell = 5 x 5 um; NVP particle radius ~4 um.  The mesh is
##   unit-agnostic (coords read 5.0); these material params set the um scale.
##   Energy release rate G_Ic in MPa*um (= uN/um).
## Catholyte / solid electrolyte = NACS (Na2Al1.35Cl4.05S).  Blocks & vars are
##   named *_NACS / *_nacs (formerly LPS/_se -- a different material).
## ============================================================================
## --- Catholyte (SE) elasto-viscoplastic properties ---
## Rate-dependent (viscoplastic) catholyte: linear elastic + static J2 yield
## (yield from Vickers hardness, Tabor sigma_y ~ H_v/3) + Norton power-law
## creep representing the viscous flow of the plasticized glassy NACS.
##
## ASSUMED temperature presets -- PLACEHOLDERS to be RECALIBRATED to the redone
## hold-segment nanoindentation (25/60/90 C). creep_coef A is tied to the
## (currently proxy) loading timescale and MUST be refit once the real
## charge/discharge loading rate is set.
##         |  25 C rigid | 60 C soft | 90 C plastic/slurry
##   E MPa |    1000     |    500    |    100
##   H_v   |     30      |     12    |     3
##   A     |    1e-7     |    1e-5   |    1e-3
##   n     |     4       |     4     |     4
## Active set = 60 C (override on the command line to sweep T / build the map).
ymod_nacs=500              # Young Modulus of SE [MPa]
Hv_nacs=12.00              # Vickers Hardness of SE [MPa]
creep_coef=1e-5          # Norton creep leading coefficient A [MPa^-n / time] (ASSUMED)
creep_exp=4              # Norton creep stress exponent n (ASSUMED; cf. LPSCl ref n=4)
pr_nacs=0.30              # Poissons Ratio of SE

ymod_nvp=90000           # Young Modulus of NVP [MPa]
pr_nvp=0.26              # Poissons Ratio of NVP

## --- Calculated Material Properties ---
ustr_nacs=${fparse Hv_nacs / 3}                             # Ultimate Strength of SE [MPa] (Coe. = 3)
ystr_nacs=${fparse ustr_nacs / 1.2}                         # Yield Strength of SE [MPa]
plstr=${fparse (ustr_nacs - ystr_nacs) / (ymod_nacs / 10)}    # Plastic Strain of SE


## --- Boundary Conditions Properties ---
sptop=10                      # Stack Pressure [MPa]
alpha_nvp=2.9927418e-4        # Thermal expansion coefficient of NVP (8.3% expansion)


## --- Output Control ---
output_times = '0.05 0.5 1.0'


[Problem]
  type = FEProblem
  solve = true
[]
[Mesh]
        [./fmg]
                type = FileMeshGenerator
                ## Mesh (AM:NACS = 50/50 area ratio)
                file = mesh_5050.msh
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
                block = 'block_NVP'
        [../]
        [./eigenstrain_yy]
                order = FIRST
                family = MONOMIAL
                block = 'block_NVP'
        [../]
        [./total_strain_xx]
                order = FIRST
                family = MONOMIAL
                block = 'block_NVP'
        [../]
        [./total_strain_yy]
                order = FIRST
                family = MONOMIAL
                block = 'block_NVP'
        [../]
[]
[Functions]
        [./temperature_load]
                type = ParsedFunction
                # NVP loading and unloading step
                expression = 'if(t<=0.5, 180*t, (-1)*(t-0.5)*180 + 90.0)'
        [../]
        [./press_ramp]
                type = PiecewiseLinear
                x = '0    0.05  1.0'
                y = '0.0  1.0   1.0'
        [../]
        [./hf_NACS]
                type = PiecewiseLinear
                x = '0 ${plstr}'
                # Hardening Function NACS
                y = '${ystr_nacs} ${ustr_nacs}'
        [../]
[]
[Modules]
        [./TensorMechanics]
         [./Master]
           [./NVP]
             strain = FINITE
                 add_variables = true
                 eigenstrain_names = eigenstrain
                 generate_output = 'stress_xx stress_yy vonmises_stress strain_xx strain_yy'
                 block = 'block_NVP'
           [../]
           [./NACS]
                 strain = FINITE
                 add_variables = true
                 generate_output = 'stress_yy stress_xx vonmises_stress plastic_strain_xx plastic_strain_yy'
                 block = 'block_NACS'
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
                block = 'block_NVP'
                rank_two_tensor = eigenstrain
                variable = eigenstrain_yy
                index_i = 1
                index_j = 1
                execute_on = 'initial timestep_end'
        [../]
        [./eigenstrain_xx]
                type = RankTwoAux
                block = 'block_NVP'
                rank_two_tensor = eigenstrain
                variable = eigenstrain_xx
                index_i = 0
                index_j = 0
                execute_on = 'initial timestep_end'
        [../]
        [./total_strain_yy]
                type = RankTwoAux
                block = 'block_NVP'
                rank_two_tensor = total_strain
                variable = total_strain_yy
                index_i = 1
                index_j = 1
                execute_on = 'initial timestep_end'
        [../]
        [./total_strain_xx]
                type = RankTwoAux
                block = 'block_NVP'
                rank_two_tensor = total_strain
                variable = total_strain_xx
                index_i = 0
                index_j = 0
                execute_on = 'initial timestep_end'
        [../]
[]
[Contact]
        [nvp_nacs]
                primary = 'block_NVP_right'
                secondary = 'block_NACS_left'
                # This penalty control the convergence issue
                # Higher Young modulus will have higher penalty parameter or vice versa
                # Tentative  number should be in the range of the young modulus
                penalty = ${ymod_nacs}
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
        #[./right_press]
                #type = Pressure
                #variable = disp_x
                #boundary = 'block_right'
                #factor   = 0.0
        #[../]
[]
[Materials]
        [./elasticity_tensor_NVP]
                type = ComputeIsotropicElasticityTensor
                block = 'block_NVP'
                youngs_modulus = ${ymod_nvp}
                poissons_ratio = ${pr_nvp}
        [../]
        [./elastic_stress_NVP]
                # NVP cathode = isotropic linear elastic (no plasticity).
                type = ComputeFiniteStrainElasticStress
                block = 'block_NVP'
        [../]
        [./thermal_expansion_strain_NVP]
                type = ComputeThermalExpansionEigenstrain
                block = 'block_NVP'
                stress_free_temperature = 0
                thermal_expansion_coeff = ${alpha_nvp}
                temperature = temp
                eigenstrain_name = eigenstrain
        [../]
        [./elasticity_tensor_NACS]
                type = ComputeIsotropicElasticityTensor
                block = 'block_NACS'
                youngs_modulus = ${ymod_nacs}
                poissons_ratio = ${pr_nacs}
        [../]
        [./isotropic_plasticity_NACS]
                type = IsotropicPlasticityStressUpdate
                block = 'block_NACS'
                yield_stress = ${ystr_nacs}
                hardening_function = hf_NACS
        [../]
        [./power_law_creep_NACS]
                # Norton viscous flow of the plasticized catholyte (rate dependence).
                # creep_rate = coefficient * vonMises^n   (activation_energy=0 ->
                # Arrhenius term = 1; T-dependence baked into coefficient above).
                type = PowerLawCreepStressUpdate
                block = 'block_NACS'
                coefficient = ${creep_coef}        # Norton A (ASSUMED; recalibrate)
                n_exponent  = ${creep_exp}         # Norton n (ASSUMED)
                activation_energy = 0              # no temperature coupling
                max_inelastic_increment = 1e-3     # bound creep increment -> substep
        [../]
        [./radial_return_stress_NACS]
                type = ComputeMultipleInelasticStress
                tangent_operator = nonlinear
                inelastic_models = 'isotropic_plasticity_NACS power_law_creep_NACS'
                block = 'block_NACS'
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
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type -mat_mumps_icntl_24'
  petsc_options_value = 'lu       mumps                       1'
  line_search = none
  nl_max_its = 99
  nl_rel_tol = 1e-6
  # abs_tol loosened 1e-8 -> 1e-7: the force-balance residual is well converged
  # by ~1e-7 (Newton then crawls linearly because the contact/creep tangent is
  # inexact), so 1e-8 only buys wasted iterations. ~2-3x fewer NL iters; gap is
  # unchanged to ~1e-5 relative. Use 1e-6 for more speed if a sweep tolerates
  # slightly looser convergence.
  nl_abs_tol = 1e-7
  l_tol = 1e-8
  start_time = 0.0
  n_startup_steps = 1
  end_time = 1
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
  file_base = rst/E${ymod_nacs}_H${Hv_nacs}_spTop${sptop}
  [exodus]
    type = Exodus
    sync_times = '${output_times}'
  []
  [csv]
    type = CSV
    sync_times = '${output_times}'
  []
[]
[Postprocessors]
  [./Gap]
    type = PointValue
    point = '0.0 3.98 0.0'
    variable = disp_y
  [../]
  [./Temp]
     type = ElementAverageValue
     variable = temp
   [../]
[]
[VectorPostprocessors]
  [disp_xy_along_arc]
    type = NodalValueSampler
    variable = 'disp_x disp_y'
    boundary = 'block_NACS_left'
    sort_by = x
    execute_on = 'TIMESTEP_END'
    outputs = csv
  []
[]