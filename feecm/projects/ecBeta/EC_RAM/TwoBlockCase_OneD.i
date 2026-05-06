
D_LiIon = 2.8e2
##IonicDiffusion
sig_LiIon=30
t_s=10e-1
delta_t= 100e-1
fieldAmp = 0.4
#n=150
factor=10
freq=${fparse 2*pi*3e5}
D_LiAtom=${fparse factor * D_LiIon}
i0_a=0.1
F = 96485
R = 8.3145
U=0
[Mesh]
    [gen]
      type = GeneratedMeshGenerator
      dim = 1
      nx = 400
      xmax = 200
      #ny = 40
      #ymax = 1
    []
    [subdomain1]
      input = gen
      type = SubdomainBoundingBoxGenerator
      bottom_left = '80 0 0'
      top_right = '200 0 0'
      block_id = 1
    []
    [interface]
      type = SideSetsBetweenSubdomainsGenerator
      input = subdomain1
      primary_block = '0'
      paired_block = '1'
      new_boundary = 'Interface'
    []
  []
  [Variables]
    [./Li_atom]
      block = 0
    [../]
    [./Li_ion]
      block = 1
    [../]
    [./pot]
      #block = 0
    [../]
  
  []
    [Functions]
      [./periodic_Pot]
        type = ParsedFunction
        value = '(fieldAmp*cos(freq*t))'
        vars = 'fieldAmp freq'
        vals = '${fieldAmp} ${freq}'  # Replace with your desired value for uValue
      [../]
        [guassian]
          type = ParsedFunction
          expression = 'amp*exp(-(x-mu)^2/(2*sig^2))'
          symbol_names = 'amp mu sig'
          symbol_values = '1e-3 150 20'
        []
      [./concDist]
        type = ParsedFunction
        expression = 'c0*(tanh(2*(x-x0))+1)'
        symbol_names = 'c0 x0'
        symbol_values = '0.225 20'  # Adjust values as needed
      [../]
      [./rect_pulse]
        type = ParsedFunction
        expression = 'if(t < start_time,0, if(t<=start_time+width,amplitude,0))'
        symbol_names = 'amplitude start_time width'
        symbol_values = '${fieldAmp} ${t_s} ${delta_t}'              # Amplitude of the pulse
      [../]
    []
    [ICs]
      [./Li_ionInitial]
         type = FunctionIC
         function = guassian
         variable = Li_ion
      [../]
      [./opIC]
        type = ConstantIC
        variable = Li_atom
        value =1e-8  
      [../]
    []
    [Kernels]
      #inactive = "DrivingForce"
      [./timeLi_ionTime]
        type = TimeDerivative
        variable = Li_ion
        block = 1
      [../]
    [./MassFluxLiIon]
        type = ADNernstPlanckConvection
        variable = Li_ion
        Voltage = pot
        diffusivity= ${D_LiIon}
        scale = 1
        block = 1
        zIons=-1
    [../]
    [./CurrentTransportLiIon]
        type = ADMigrationDiffusion
        variable = pot
        c = Li_ion
        #T= T
        c0=1e-4
        conductivity = ${sig_LiIon}
        scale = 1
        block = 1
        zIons = -1
    [../]
      [./timeLi_AtomTime]
          type = TimeDerivative
          variable = Li_atom
          block = 0
      [../]
      [./DiffusionLiAtom]
          type = PhaseFieldLaplace
          variable = Li_atom
          k0= ${fparse D_LiAtom}
          scale = 1
          block = 0
      [../]
  
    []
    [AuxVariables]
      [./E_field]
        order = CONSTANT
        family = MONOMIAL
      [../]
      [totalLiConc][]
    []
    
    [AuxKernels]
      [./E_fieldAux]
        type = FunctionAux
        variable = E_field
        function = rect_pulse
      [../]
      [u]
          type = ParsedAux
          variable = totalLiConc
          coupled_variables = 'Li_atom'
          expression = 'Li_atom'
          block = 0
      []
      [v]
          type = ParsedAux
          variable = totalLiConc
          coupled_variables = 'Li_ion'
          expression = 'Li_ion'
          block = 1
      []
    []
    [Materials]
        [charge_transfer_anode_elyte]
            type = ChargeTransferReaction
            #electrode = true
            charge_transfer_current_density = ibv
            charge_transfer_mass_flux = je
            charge_transfer_heat_flux = he
            electric_potential = pot
            neighbor_electric_potential = pot
            charge_transfer_coefficient = 0.5
            exchange_current_density = ${i0_a}
            faraday_constant = ${F}
            ideal_gas_constant = ${R}
            temperature = 300
            open_circuit_potential = ${U}
            boundary = 'Interface'
          []
    []
    [InterfaceKernels]
      active = 'diffusion'
      [./diffusion]
        type = InterfaceDiffusion
        variable = Li_atom
        neighbor_var = Li_ion
        boundary = Interface
        D = ${fparse D_LiAtom}
        D_neighbor = ${fparse D_LiIon}
      [../]
    [./penalty]
        type = PenaltyInterfaceDiffusion
        variable = Li_atom
        neighbor_var = Li_ion
        boundary = Interface
        penalty = 1e3
    [../]
    # [negative_current]
    #     type = MaterialInterfaceNeumannBC
    #     variable = pot
    #     neighbor_var = pot
    #     prop = ibv
    #     factor = 1
    #     boundary = 'Interface'
    # []
    []
    [Postprocessors]
      [potSE]
        type = ElementAverageValue
        variable = pot
      []
      [appliedField]
        type = ElementAverageValue
        variable = E_field
      []
      [Conc_LiIon]
        type = ElementAverageValue
        variable = Li_ion
        block = 1
      []
      [Conc_LiAtom]
          type = ElementAverageValue
          variable = Li_atom
          block = 0
      []
    []
    [Preconditioning]
      [./smp]
        type = SMP
        full = true
        petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -snes_atol -snes_rtol -ksp_atol -ksp_rtol'
        petsc_options_value = 'hypre boomeramg 100 1e-8 1e-8 1e-8 1e-8'
      [../]
    []
    
    [Executioner]
      type = Transient
      solve_type = NEWTON
      automatic_scaling = true
      start_time = 0.0
      #dtmin = 1e-10
      # dtmax=5e-8
      end_time=50
      l_max_its = 50
      nl_max_its = 20
      nl_rel_tol = 1e-6
      nl_abs_tol = 1e-8
      [./TimeStepper]
        type = IterationAdaptiveDT
        optimal_iterations = 6
        iteration_window = 2
        growth_factor = 1.5
        cutback_factor = 0.5
        dt = 1e-4
      [../]
      # [./TimeStepper]
      #     type = ConstantDT
      #     dt = 1e-1
      #     [../]
          dtmin = 1e-4
    []
    [BCs]
      [ConcLiFluxIon]
          type=NeumannBC
          variable=Li_ion
          boundary='right'
          value ='0.0'
      []
      [ConcLiIon]
        type=DirichletBC
        variable=Li_ion
        boundary='right'
        value ='0.0'
      []
      [ConcLiLeft]
          type=NeumannBC
          variable=Li_atom
          boundary='left'
          value ='0'
      []
      # [Flux]
      #   type=VacuumBC
      #   variable=conc
      #   boundary='left right top bottom'
      #   # value ='0'
      # []
      # [ConcLeft]
      #   type=NeumannBC
      #   variable=conc
      #   boundary='left right'
      #   value ='0.0'
      # []
    #   [LeftPot]
    #     type=DirichletBC
    #     variable=pot
    #     boundary='Interface'
    #     value ='0.0'
    #   []
      # [top_bc]
      #   type = ButlerVolmerIonics
      #   variable = pot
      #   boundary = Interface
      #   LiPotRef = 0
      #   #LiPotRef = -2
      #   ex_current= 0.5
      # []
      [InterfacePot]
          type=DirichletBC
          variable=pot
          boundary='Interface'
          value ='0'
        []
        [RightPot]
          type=DirichletBC
          variable=pot
          boundary='right'
          value = ${fieldAmp}
          #function =rect_pulse
        []
    
    []
    
    [Outputs]
      [./out]
        type = Exodus
        file_base = TwoBlockCase_bvOFF_${fieldAmp}  #twoBlockCase_${fieldAmp}_NeuBc_EFturnOn
        elemental_as_nodal = true
      [../]
    []
    