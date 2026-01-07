/* ----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine
    http://www.dcs-computing.com, office@dcs-computing.com
    LIGGGHTS® is part of CFDEM®project: http://www.liggghts.com | http://www.cfdem.com

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    Created: Joseph Vazquez Mercado, RIT 2025
    Copyright 2024-     DCS Computing GmbH, Linz
    Notes: Radius expansion header for metal-carbon lithium alloying system
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(metal_radius_expansion,FixMetalRadiusExpansion)

#else

#ifndef LMP_FIX_METAL_RADIUS_EXPANSION_H
#define LMP_FIX_METAL_RADIUS_EXPANSION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMetalRadiusExpansion : public Fix {
 public:
  FixMetalRadiusExpansion(class LAMMPS *, int, char **);
  virtual ~FixMetalRadiusExpansion();
  virtual int setmask();
  virtual void init();
  virtual void setup(int);
  virtual void post_create();
  virtual void end_of_step();  
  
  double compute_scalar();
  double compute_vector(int);
  void trigger_radius_update();
  
 protected:
  // Volume parameters for metal-Li alloying
  double Omega_metal;        // Molar volume of metal (10.96e-6 m³/mol for Ag)
  double Omega_Li_eff;       // Effective molar volume of Li in alloy (9.0e-6 m³/mol)
  
  // Update control (no max limit)
  int update_frequency;      // Steps between radius updates
  int update_count;          // Current update count (for diagnostics)
  int next_update_step;      // Next step to update radius
  bool manual_trigger;       // Manual trigger for radius update
  
  // Property pointers
  double *lithium_content;   // Li/Metal molar ratio
  double *li_mols;           // Moles of lithium in particle
  double *metal_mols;        // Moles of metal in particle (constant)
  double *initial_radius;    // Initial particle radius (μm)
  double *particle_volume;   // Current particle volume (m³)
  
  // Fix pointers
  class FixPropertyAtom *fix_lithium_content;
  class FixPropertyAtom *fix_li_mols;
  class FixPropertyAtom *fix_metal_mols;
  class FixPropertyAtom *fix_initial_radius;
  class FixPropertyAtom *fix_particle_volume;
  
  // Particle types
  int metal_type;   // Atom type for metal particles (default: 1)
  int carbon_type;  // Atom type for carbon particles (default: 2)
  
  // Methods
  void update_particle_radius();
  void updatePtrs();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix metal_radius_expansion command

Self-explanatory. Check the input script syntax and compare to the
documentation for the command.

E: Fix metal_radius_expansion requires lithiumContent property

The lithiumContent property/atom must be available.

E: Fix metal_radius_expansion requires lithiumMols property

The lithiumMols property/atom must be available (created by fix metal_carbon_diffusion).

E: Fix metal_radius_expansion requires metalMols property

The metalMols property/atom must be available (created by fix metal_carbon_diffusion).

E: Could not find required property/atom fixes

Required property fixes (initialRadius, particleVolume) were not found.

*/