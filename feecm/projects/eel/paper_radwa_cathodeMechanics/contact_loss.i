## Young Modulus of LPSCL MPa
ymod=1000
## Yield Strength of LPSCL MPa
ystr=196
ystr_plus_4=${fparse ystr+4.0}
## Stack Pressure MPa
sptop=6
[Problem]
  type = FEProblem
  solve = true
[] 
[Mesh] 
        [./fmg] 
                type = FileMeshGenerator 
                ## Mesh file with NMC:LPS=70/30 ratio
                file = data.msh 
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
                y = '12600'
        [../] 
        [./hf_LPS] 
                type = PiecewiseLinear 
                x = '0 0.004'
                # Hardening Function LPS
                y = '${ystr} ${ystr_plus_4}' 
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
                penalty = 1e6 #100
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
                youngs_modulus = 177500 
                poissons_ratio = 0.33 
        [../] 
        [./isotropic_plasticity_NMC] 
                type = IsotropicPlasticityStressUpdate 
                block = 'block_NMC' 
                yield_stress = 12600.0 
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
                ## use 1/3 of 8.1 volumetric expansion which would give the following thermal expansion coeff

                thermal_expansion_coeff = 0.00028     #0.00075 
                temperature = temp 
                eigenstrain_name = eigenstrain 
        [../] 
        [./elasticity_tensor_LPS] 
                type = ComputeIsotropicElasticityTensor 
                block = 'block_LPS' 
                youngs_modulus = ${ymod} #25300.0  
                poissons_ratio = 0.23
        [../] 
        [./isotropic_plasticity_LPS] 
                type = IsotropicPlasticityStressUpdate 
                block = 'block_LPS' 
                yield_stress = ${ystr} #5000.0 
                hardening_function = hf_LPS
                #hardening_constant = 50 
                # relative_tolerance = 1e-6
                # absolute_tolerance = 1e-6
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
        nl_max_its = 99
        nl_rel_tol = 1e-7 #1e-7 
        nl_abs_tol = 1e-9 #1e-9 
        l_tol = 1e-8 #1e-8
        start_time = 0.0 
        n_startup_steps = 1 
        end_time = 1 
        # dt = 0.025 
        # dtmin = 1e-3 #0.001 
        [TimeStepper]
            type = IterationAdaptiveDT
            optimal_iterations = 1
            linear_iteration_ratio = 1
            dt = 2e-2
        []
[] 
[Outputs]
  exodus = true
  file_base = rst/ystr${ystr}_E${ymod}_spTop${sptop}
  [./csv]
    type = CSV
    execute_on = 'final'
  [../]
[]
[Postprocessors]
  [./Gap]
    type = PointValue
    point = '0.0 0.00250001 0.0'
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
