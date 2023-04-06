!config navigation breadcrumbs=False scrollspy=False

# Electro_chemo_mechanics App

The Electro-Chemo-Mechanics (ECM) App is a MOOSE based application, with a set of tools that 
solves the continuum mechanics problem of Li deposition on a solid state electrolyte, typically
used in an All Solid State Battery (ASSB). It provides a simple approach with a plugable 
architecture that can be used to implement advanced material and interface models 

The key characteristics of the ECM App are 

- Li Plasticity (ref)
- Large deformation interface mechanics 

    - Using mortar FEM and electro-chemistry (Butler-Volmer Kinetics) 
    
    - Using Cohesive Zone model (CZM) with electro-chemistry 

