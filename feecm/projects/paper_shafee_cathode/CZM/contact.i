## --- Bulk Material Properties ---
ymod_se=100             # Young Modulus of SE [MPa] (Sweeping from 100 to 1000 MPa)
Hv_se=3.00              # Vickers Hardness of SE [MPa] (Sweeping from 3 to 30 MPa)
pr_se=0.30              # Poissons Ratio of SE

ustr_am=10000           # Ultimate Strength of AM [MPa] (Some High Value to Ignore AM Plasticity)
ymod_am=90000           # Young Modulus of AM [MPa]
pr_am=0.26              # Poissons Ratio of AM


## --- Calculated Material Properties ---
ustr_se=${fparse Hv_se / 3}                             # Ultimate Strength of SE [MPa] (Coe. = 3)
ystr_se=${fparse ustr_se / 1.2}                         # Yield Strength of SE [MPa]
plstr=${fparse (ustr_se - ystr_se) / (ymod_se / 10)}    # Plastic Strain of SE


## --- Boundary Conditions Properties ---
sptop=10                      # Stack Pressure [MPa]
alpha_nmc=2.9927418e-4        # Thermal expansion coefficient of NMC (8.3% expansion)


[Problem]
  type = FEProblem
  solve = true
[]
[Mesh]
        [./fmg]
                type = FileMeshGenerator
                ## Mesh file with AM:SE=50/50 volume ratio
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
                # NMC loading and unloading step
                expression = 'if(t<=0.5, 180*t, (-1)*(t-0.5)*180 + 90.0)'
        [../]
        [./hf_NMC]
                type = PiecewiseLinear
                x = '0'
                y = '${ustr_am}'
        [../]
        [./hf_LPS]
                type = PiecewiseLinear
                x = '0 ${plstr}'
                # Hardening Function LPS
                y = '${ystr_se} ${ustr_se}'
        [../]
[]
[Modules]
        [./TensorMechanics]
         [./Master]
           [./NMC]
             strain = FINITE
                 add_variables = true
                 eigenstrain_names = eigenstrain
                 generate_output = 'stress_xx stress_yy vonmises_stress strain_xx strain_yy'
                 block = 'block_NMC'
           [../]
           [./LPS]
                 strain = FINITE
                 add_variables = true
                 generate_output = 'stress_yy stress_xx vonmises_stress plastic_strain_xx plastic_strain_yy'
                 block = 'block_LPS'
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
        [nmc_lps]
                primary = 'block_NMC_right'
                secondary = 'block_LPS_left'
                # This penalty control the convergence issue
                # Higher Young modulus will have higher penalty parameter or vice versa
                # Tentative  number should be in the range of the young modulus
                penalty = ${ymod_se}
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
        [../]
        #[./right_press]
                #type = Pressure
                #variable = disp_x
                #boundary = 'block_right'
                #factor   = 0.0
        #[../]
[]
[Materials]
        [./elasticity_tensor_NMC]
                type = ComputeIsotropicElasticityTensor
                block = 'block_NMC'
                youngs_modulus = ${ymod_am}
                poissons_ratio = ${pr_am}
        [../]
        [./isotropic_plasticity_NMC]
                type = IsotropicPlasticityStressUpdate
                block = 'block_NMC'
                yield_stress = ${ustr_am}
                hardening_function = hf_NMC
        [../]
        [./radial_return_stress_NMC]
                type = ComputeMultipleInelasticStress
                tangent_operator = elastic
                inelastic_models = 'isotropic_plasticity_NMC'
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
        [./elasticity_tensor_LPS]
                type = ComputeIsotropicElasticityTensor
                block = 'block_LPS'
                youngs_modulus = ${ymod_se}
                poissons_ratio = ${pr_se}
        [../]
        [./isotropic_plasticity_LPS]
                type = IsotropicPlasticityStressUpdate
                block = 'block_LPS'
                yield_stress = ${ystr_se}
                hardening_function = hf_LPS
        [../]
        [./radial_return_stress_LPS]
                type = ComputeMultipleInelasticStress
                tangent_operator = elastic
                inelastic_models = 'isotropic_plasticity_LPS'
                block = 'block_LPS'
        [../]
[]
[Executioner]
        type = Transient
        automatic_scaling = true
        solve_type = NEWTON
        petsc_options_iname = -pc_type
        petsc_options_value = lu
        line_search = none
        nl_max_its = 99
        nl_rel_tol = 1e-6
        nl_abs_tol = 1e-9
        l_tol = 1e-8
        start_time = 0.0
        n_startup_steps = 1
        end_time = 1
        [TimeStepper]
            type = IterationAdaptiveDT
            dt = 0.01
            optimal_iterations = 10
            growth_factor = 1.5
            cutback_factor = 0.5
            cutback_factor_at_failure = 0.2
            linear_iteration_ratio = 100
        []
[]
[Outputs]
  exodus = true
  file_base = rst/ystr${ystr_se}_E${ymod_se}_spTop${sptop}
  [./csv]
    type = CSV
    execute_on = 'final'
  [../]
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
[./disp_xy_along_arc]
    type = NodalValueSampler
    variable = 'disp_x disp_y'
    boundary = 'block_LPS_left'
    sort_by = x
    execute_on = 'final'
  [../]
[]