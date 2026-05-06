/* ----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine
    Contributing author: Joseph Vazquez Mercado, RIT 2025
    Copyright 2024-     DCS Computing GmbH, Linz
    
    Implicit (Backward Euler) Lithium Diffusion with Chemo-Mechanical Coupling
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
FixStyle(lithium_diffusion_implicit,FixLithiumDiffusionImplicit)
#else

#ifndef LMP_FIX_LITHIUM_DIFFUSION_IMPLICIT_H
#define LMP_FIX_LITHIUM_DIFFUSION_IMPLICIT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPropertyAtom;
class FixPropertyAtomLithiumContent;
class NeighList;

class FixLithiumDiffusionImplicit : public Fix {
 public:
  FixLithiumDiffusionImplicit(class LAMMPS *, int, char **);
  ~FixLithiumDiffusionImplicit();

  void post_create();
  int setmask();
  void init();
  void setup(int);
  void init_list(int, NeighList *);
  void end_of_step();
  double compute_scalar();
  double compute_vector(int);

 private:
  void updatePtrs();
  void ensure_arrays_sized();
  void implicit_diffusion_step(double dt_diff);
  double calculate_contact_area(int, int);
  double get_diffusion_coefficient(int, int);
  double get_c_li_max(int);
  double get_V_exp_max(int);
  double get_Omega_Li(int);
  double compute_mobility(double, double, double);
  double compute_effective_potential(int);
  double compute_D_CB(double x_norm);
  double compute_D_particle(int i);
  void   update_per_atom_D();

  // Physical constants
  double R, T, F;

  // Max concentrations (mol/m³)
  double c_li_max_LM, c_li_max_AM, c_li_max_CB;

  // Volume expansion
  double V_exp_max_AM, V_exp_max_CB, V_exp_max_LM;

  // Partial molar volumes of Li in each phase (m³/mol)
  double Omega_Li_AM;          // Active material
  double Omega_Li_CB;          // Carbon black
  double Omega_Li_LM;          // Lithium metal

  // Lithium content bounds
  double initial_lithium_content, target_lithium_content, max_lithium_content;

  // Diffusion coefficients
  //   AM, LM: constant (user-specified)
  //   CB:     SOC-dependent from GITT rational fit
  double D_AM;              // m²/s — constant for active material
  double D_LM;              // m²/s — constant for lithium metal
  double D_CB_floor;        // m²/s — safety floor for CB rational fit

  // Implicit solver parameters
  double diffusion_dt;
  int    update_every;
  int    max_iter;
  double tol;
  int    num_substeps;

  // Diagnostics
  int    last_iter_count;
  double last_residual;
  double total_diffusion_time;

  // Temporary arrays for Jacobi iteration
  int    nmax_alloc;
  double *c_old;
  double *c_new;
  double *cross_flux_arr;

  // Per-atom property pointers
  double *lithium_content, *lithium_concentration, *equilibrium_potential;
  double *hydrostatic_stress;
  double *current_el, *diffusion_coefficient, *lithium_flux;
  double *li_mols, *initial_volume;

  // Fix pointers
  FixPropertyAtom *fix_initial_volume, *fix_lithium_content;
  FixPropertyAtom *fix_lithium_concentration, *fix_equilibrium_potential;
  FixPropertyAtom *fix_hydrostatic_stress;
  FixPropertyAtom *fix_current_el, *fix_diffusion_coefficient;
  FixPropertyAtom *fix_lithium_flux, *fix_li_mols;
  FixPropertyAtomLithiumContent *fix_lithium_content_manager;

  int AM_type, CB_type, LM_type;
  NeighList *list;
};

}

#endif
#endif