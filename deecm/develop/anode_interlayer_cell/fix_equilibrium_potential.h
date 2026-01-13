/* ----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine
    Contributing author: Joseph Vazquez Mercado, RIT 2025
    Copyright 2024-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
FixStyle(equilibrium_potential,FixEquilibriumPotential)
#else

#ifndef LMP_FIX_EQUILIBRIUM_POTENTIAL_H
#define LMP_FIX_EQUILIBRIUM_POTENTIAL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPropertyAtom;

class FixEquilibriumPotential : public Fix {
 public:
  FixEquilibriumPotential(class LAMMPS *, int, char **);
  ~FixEquilibriumPotential();
  
  void post_create();
  int setmask();
  void init();
  void setup(int);
  void pre_force(int);
  void post_force(int);
  double compute_scalar();
  
  // Material type enumeration
  enum { MATERIAL_AG, MATERIAL_C, MATERIAL_NMC811 };

 private:
  void updatePtrs();
  void calculate_equilibrium_potential();
  double compute_OCV(double x_li, double x_max, int material);
  
  // Physical constants
  double R;              // Gas constant (J/(mol·K))
  double T;              // Temperature (K)
  double F;              // Faraday constant (C/mol)
  
  // Material parameters
  double x_Li_max_AM;    // Max Li content for AM
  double x_Li_max_CB;    // Max Li content for CB
  double U_OCV_LM;       // OCV for lithium metal (V vs Li/Li+)
  
  // Material types for AM and CB
  int AM_material;       // Material formula to use for AM particles
  int CB_material;       // Material formula to use for CB particles
  
  // Particle type identifiers
  int AM_type;
  int CB_type;
  int LM_type;
  
  // Per-atom properties
  double *lithium_content;
  double *equilibrium_potential;
  
  // Fix pointers
  FixPropertyAtom *fix_lithium_content;
  FixPropertyAtom *fix_equilibrium_potential;
};

}

#endif
#endif