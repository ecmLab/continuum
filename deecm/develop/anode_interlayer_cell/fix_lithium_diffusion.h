/* ----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine
    Contributing author: Joseph Vazquez Mercado, RIT 2025
    Copyright 2024-     DCS Computing GmbH, Linz
    Notes: Multi-material Li diffusion: AM, CB, and LM (constant source)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(lithium_diffusion,FixLithiumDiffusion)

#else

#ifndef LMP_FIX_LITHIUM_DIFFUSION_H
#define LMP_FIX_LITHIUM_DIFFUSION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPropertyAtom;
class FixPropertyAtomLithiumContent;
class NeighList;

class FixLithiumDiffusion : public Fix {
 public:
  FixLithiumDiffusion(class LAMMPS *, int, char **);
  virtual ~FixLithiumDiffusion();
  void post_create();
  int setmask();
  void init();
  void setup(int);
  void init_list(int, class NeighList *);
  void pre_force(int);
  void post_force(int);
  double compute_scalar();
  void updatePtrs();

 protected:
  // Physical constants
  double F;                      // Faraday constant (96485 C/mol)
  // Material-specific maximum concentrations (mol/m^3)
  double c_li_max_LM;
  double c_li_max_AM;
  double c_li_max_CB;

  double V_exp_max_AM;          // Max volume expansion for AM
  double V_exp_max_CB;          // Max volume expansion for CB
  double V_exp_max_LM;          // Max volume expansion for LM
  
  // Molar volumes (m³/mol)
  double Omega_Li_LM;            // Li molar volume in LM (12.97e-6 m³/mol)

  // Lithium content parameters
  double initial_lithium_content;
  double target_lithium_content;
  double max_lithium_content;

  // Diffusion coefficients (m²/s) for each material pair
  double D_AM_AM;                // AM to AM diffusion
  double D_AM_CB;                // AM to CB diffusion
  double D_AM_LM;                // AM to LM diffusion
  double D_CB_AM;                // CB to AM diffusion
  double D_CB_CB;                // CB to CB diffusion
  double D_CB_LM;                // CB to LM diffusion

  // Property pointers
  double *lithium_content;
  double *lithium_concentration;
  double *current_SE_Li;
  double *diffusion_coefficient;
  double *lithium_flux;
  double *li_mols;
  double *initial_volume;              // Store initial particle volumes
  
  // Fix pointers
  class FixPropertyAtom *fix_lithium_content;
  class FixPropertyAtom *fix_lithium_concentration;
  class FixPropertyAtom *fix_current_SE_Li;
  class FixPropertyAtom *fix_diffusion_coefficient;
  class FixPropertyAtom *fix_lithium_flux;
  class FixPropertyAtom *fix_li_mols;
  class FixPropertyAtomLithiumContent *fix_lithium_content_manager;
  class FixPropertyAtom *fix_initial_volume; // Fix pointer for the property
  
  // Particle types
  int AM_type;                   // Active Material type (default 1)
  int CB_type;                   // Carbon Black type (default 2)
  int LM_type;                   // Lithium Metal type (default 4)
  
  // Neighbor list
  class NeighList *list;
  
  // Methods
  void update_lithium_content();
  double calculate_contact_area(int, int);
  double get_diffusion_coefficient(int, int);
  double get_c_li_max(int); // New helper function
  double get_V_exp_max(int); // New helper function
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix lithium_diffusion command

Self-explanatory. Check the input script syntax and compare to the
documentation for the command.

E: Fix lithium_diffusion requires lithiumContent property

The lithiumContent property/atom must be defined before fix lithium_diffusion.

E: Fix lithium_diffusion requires lithiumConcentration property

The lithiumConcentration property/atom must be defined before fix lithium_diffusion.

E: Could not find required property/atom fixes

One or more of the required property/atom fixes (diffusionCoefficient,
lithiumFlux, lithiumMols) could not be found or created.

*/