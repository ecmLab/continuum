/* ----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine
    Contributing author: Joseph Vazquez Mercado, RIT 2025
    Copyright 2024-     DCS Computing GmbH, Linz
    Notes: Radius expansion for AM and CB particles based on Li content
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(radius_expansion,FixRadiusExpansion)

#else

#ifndef LMP_FIX_RADIUS_EXPANSION_H
#define LMP_FIX_RADIUS_EXPANSION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRadiusExpansion : public Fix {
 public:
  FixRadiusExpansion(class LAMMPS *, int, char **);
  virtual ~FixRadiusExpansion();
  virtual int setmask();
  virtual void init();
  virtual void setup(int);
  virtual void post_create();
  virtual void end_of_step();  
  
  double compute_scalar();
  double compute_vector(int);
  void trigger_radius_update();
  
 protected:
  // Volume parameters - effective molar volumes (m³/mol)
  double Omega_Li_AM;        // Effective molar volume of Li in AM
  double Omega_Li_CB;        // Effective molar volume of Li in CB
  
  // Update control
  int update_frequency;      // Steps between radius updates
  int max_updates;           // Maximum number of radius updates
  int update_count;          // Current update count
  int next_update_step;      // Next step to update radius
  bool manual_trigger;       // Manual trigger for radius update
  
  // Property pointers
  double *lithium_content;
  double *li_mols;           // Moles of lithium (from lithium_diffusion)
  double *particle_volume;
  double *initial_volume;   // Initial particle volume
  
  // Fix pointers
  class FixPropertyAtom *fix_lithium_content;
  class FixPropertyAtom *fix_li_mols;
  class FixPropertyAtom *fix_particle_volume;
  class FixPropertyAtom *fix_initial_volume;
  
  // Particle types
  int AM_type;               // Active Material type (default 1)
  int CB_type;               // Carbon Black type (default 2)
  
  // Methods
  void update_particle_radius();
  void updatePtrs();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix radius_expansion command

Self-explanatory. Check the input script syntax and compare to the
documentation for the command.

E: Fix radius_expansion requires lithiumContent property

The lithiumContent property/atom must be defined before fix radius_expansion.

E: Fix radius_expansion requires lithiumMols property

The lithiumMols property/atom must be created by fix lithium_diffusion
before fix radius_expansion is defined.

E: Could not find required property/atom fixes

One or more required property/atom fixes (initialRadius, particleVolume)
could not be found or created.

*/