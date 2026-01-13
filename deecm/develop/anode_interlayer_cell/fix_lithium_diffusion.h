/* ----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine
    Contributing author: Joseph Vazquez Mercado, RIT 2025
    Copyright 2024-     DCS Computing GmbH, Linz
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
  ~FixLithiumDiffusion();
  
  void post_create();
  int setmask();
  void init();
  void setup(int);
  void init_list(int, NeighList *);
  void pre_force(int);
  void post_force(int);
  double compute_scalar();

 private:
  void updatePtrs();
  void update_lithium_content();
  double calculate_contact_area(int, int);
  double get_diffusion_coefficient(int, int);
  double get_c_li_max(int);
  double get_V_exp_max(int);
  double compute_mobility(double c_i, double c_j, double D_ij);
  
  // Physical constants
  double R;                      // Gas constant (J/(mol·K))
  double T;                      // Temperature (K)
  double F;                      // Faraday constant (C/mol)
  
  // Maximum lithium concentrations (mol/m³)
  double c_li_max_LM;
  double c_li_max_AM;
  double c_li_max_CB;
  
  // Volume expansion parameters
  double V_exp_max_AM;
  double V_exp_max_CB;
  double V_exp_max_LM;
  
  // Molar volume
  double Omega_Li_LM;
  
  // Lithium content parameters
  double initial_lithium_content;
  double target_lithium_content;
  double max_lithium_content;
  
  // Diffusion coefficients (m²/s)
  // Same-type: concentration-driven
  double D_AM_AM;
  double D_CB_CB;
  // Cross-type: used for mobility calculation
  double D_AM_CB;
  double D_AM_LM;
  double D_CB_AM;
  double D_CB_LM;
  
  // Per-atom property pointers
  double *lithium_content;
  double *lithium_concentration;
  double *equilibrium_potential;
  double *current_SE_Li;
  double *diffusion_coefficient;
  double *lithium_flux;
  double *li_mols;
  double *initial_volume;
  
  // Fix pointers
  FixPropertyAtom *fix_initial_volume;
  FixPropertyAtom *fix_lithium_content;
  FixPropertyAtom *fix_lithium_concentration;
  FixPropertyAtom *fix_equilibrium_potential;
  FixPropertyAtom *fix_current_SE_Li;
  FixPropertyAtom *fix_diffusion_coefficient;
  FixPropertyAtom *fix_lithium_flux;
  FixPropertyAtom *fix_li_mols;
  FixPropertyAtomLithiumContent *fix_lithium_content_manager;
  
  // Particle type identifiers
  int AM_type;
  int CB_type;
  int LM_type;
  
  // Neighbor list
  NeighList *list;
};

}

#endif
#endif