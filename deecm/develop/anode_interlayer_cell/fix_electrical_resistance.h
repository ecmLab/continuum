/* ----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine
    Contributing author: Based on Zhang et al. (2024) Powder Technology
    "An electro-mechanical contact model for particulate systems"
    
    Implements coupled bulk + contact resistance model for battery
    particle systems with full MPI support.
    
    Uses direct geometric calculations from LIGGGHTS particle data.
    All internal calculations performed in SI units.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
FixStyle(electrical_resistance, FixElectricalResistance)
#else

#ifndef LMP_FIX_ELECTRICAL_RESISTANCE_H
#define LMP_FIX_ELECTRICAL_RESISTANCE_H

#include "fix.h"
#include <vector>

namespace LAMMPS_NS {

class FixPropertyAtom;
class PairGran;
class NeighList;

class FixElectricalResistance : public Fix {
 public:
  FixElectricalResistance(class LAMMPS *, int, char **);
  ~FixElectricalResistance();

  void post_create();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void pre_force(int);
  void post_force(int);
  double compute_scalar();
  double compute_vector(int);

  // Communication functions for ghost atoms
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);

 private:
  // Helper functions
  void updatePtrs();
  void calculate_resistances();
  
  // Direct geometric calculations (all in SI units)
  double compute_overlap(int i, int j, double &L_ij);
  double compute_contact_radius(double overlap, double r_i, double r_j);
  double compute_contact_area(double overlap, double r_i, double r_j);
  
  // Resistance calculations (all in SI units)
  double compute_contact_resistance(double rho_i, double rho_j, double r_c);
  double compute_bulk_resistance(double rho, double r_i, double r_j, double L_ij);
  double compute_total_resistance(int i, int j);

  // Property access
  double get_resistivity(int type);
  double get_effective_radius(double r_i, double r_j);

  // Physical constants
  static constexpr double PI = 3.14159265358979323846;
  
  // Unit conversion: LIGGGHTS micro units to SI
  // LIGGGHTS uses μm for length in micro unit system
  static constexpr double LENGTH_CONV = 1.0e-6;  // μm to m
  
  // Temperature
  double T_;                    // Temperature [K]

  // Particle type identifiers
  int AM_type_;                 // Active material type (Metal)
  int CB_type_;                 // Carbon black type
  int SE_type_;                 // Solid electrolyte type (boundary)
  int LM_type_;                 // Lithium metal type (boundary)
  int wall_type_;               // Wall type
  int ntypes_;                  // Number of atom types

  // Material properties - stored as resistivity [Ω·m] (SI units)
  // Input is conductivity [S/m], converted: ρ = 1/σ
  std::vector<double> resistivity_;
  double wall_resistivity_;

  // Minimum overlap threshold [m] (converted from input in μm)
  double min_overlap_;

  // Per-atom property pointers
  double *total_conductance_;       // Sum of conductances to all neighbors [S]
  double *total_resistance_;        // Effective resistance (1/conductance) [Ω]
  double *contact_count_;           // Number of contacts per particle
  double *avg_contact_resistance_;  // Average contact resistance [Ω]

  // Fix pointers for per-atom properties
  FixPropertyAtom *fix_total_conductance_;
  FixPropertyAtom *fix_total_resistance_;
  FixPropertyAtom *fix_contact_count_;
  FixPropertyAtom *fix_avg_contact_resistance_;

  // Pair style and neighbor list
  PairGran *pair_gran_;
  NeighList *list_;
  int irequest_;

  // Newton pair flag
  int newton_pair_;

  // Communication sizes
  int nmax_comm_;
  
  // Output values (all in SI: Ω for resistance, S for conductance)
  double total_system_resistance_;
  double avg_pair_resistance_;
  int total_contact_pairs_;
};

}

#endif
#endif