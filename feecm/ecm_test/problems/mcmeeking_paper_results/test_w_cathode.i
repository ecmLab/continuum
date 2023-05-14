# Test to check results of mcmeeking paper
# Test overpotential on the surface etc.

current_density = 10e-3
k_anode = 1e-3
faraday = 96.4853329
temperature = 298
gas_constant = 8.31446
c_ref = 0
c_ref_pore = 2.5925e-2
k_pore = 1e-3
c_init_Li = 0
c_init_C = 2.5925e-2
c_max_C = 3.05e-2
pressure = 1e-4
end_time = 150
[Mesh]
  construct_node_list_from_side_list = true
  [./mesh]
    type = FileMeshGenerator
    file = 'data/test_w_cathode.e'
  [../]
  # ------------ Interfaces ---------------
  # --- interface between interlayer and ceramic
  [./secondary_interLayer_block]
    type = LowerDBlockFromSidesetGenerator
    input = mesh
    sidesets = 'interLayer_bottom'
    new_block_name = 'interLayer_s_block'
  [../]
  [./primary_interLayer_block]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_interLayer_block
    sidesets = 'blockCeramic_top'
    new_block_name = 'interLayer_p_block'
  [../]
  # --- interface between pore and ceramic
  [./secondary_Pore_block]
    type = LowerDBlockFromSidesetGenerator
    input = primary_interLayer_block
    sidesets = 'blockPore_bottom'
    new_block_name = 'blockPore_s_block'
  [../]
  [./primary_Pore_block]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_Pore_block
    sidesets = 'blockCeramic_top'
    new_block_name = 'blockPore_p_block'
  [../]
  # --- interface between ceramic and cathode
  [./secondary_Cathode_block]
    type = LowerDBlockFromSidesetGenerator
    input = primary_Pore_block
    sidesets = 'blockCeramic_bottom'
    new_block_name = 'blockCathode_s_block'
  [../]
  [./primary_Cathode_block]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_Cathode_block
    sidesets = 'blockCathode_top'
    new_block_name = 'blockCathode_p_block'
  [../]

[]

[GlobalParams]
  displacements = 'ux uy'
[]
[Modules/TensorMechanics/Master]
  [./all]
    add_variables = true
    strain = FINITE
    volumetric_locking_correction = true
    generate_output = 'stress_xx stress_yy strain_xx strain_yy vonmises_stress hydrostatic_stress'
    use_automatic_differentiation = true
    block = 'blockCeramic interLayer blockPore blockMetal blockCathode '
    extra_vector_tags = 'ref'
  [../]
[]

[Contact]
  [./metal_interlayer]
    formulation = mortar
    model = frictionless
    primary = 'blockMetal_bottom'
    secondary ='interLayer_top'
    normal_smoothing_distance = 0.5
    tangential_tolerance = 0.25
    mesh = primary_Cathode_block
  [../]
  [./metal_pore]
    formulation = mortar
    model = frictionless
    primary = 'blockMetal_right'
    secondary ='blockPore_left'
    normal_smoothing_distance = 0.5
    tangential_tolerance = 0.25
    mesh = metal_interlayer_secondary_subdomain_generator
  [../]
  [./interlayer_pore]
    formulation = mortar
    model = frictionless
    primary = 'interLayer_right'
    secondary ='blockPore_left'
    normal_smoothing_distance = 0.5
    tangential_tolerance = 0.25
    mesh = metal_pore_secondary_subdomain_generator
  [../]
[]

[Variables]
  [./ux]
    block = 'blockCeramic interLayer blockPore blockMetal blockCathode '
  [../]
  [./uy]
    block = 'blockCeramic interLayer blockPore blockMetal blockCathode '
  [../]
