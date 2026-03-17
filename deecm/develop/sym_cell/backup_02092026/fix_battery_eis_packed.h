/* ----------------------------------------------------------------------
    LIGGGHTS® DEM simulation engine
    
    Header for fix_battery_eis with corrected current conservation
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
FixStyle(battery/eis,FixBatteryEIS)
#else

#ifndef LMP_FIX_BATTERY_EIS_H
#define LMP_FIX_BATTERY_EIS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBatteryEIS : public Fix {
 public:
  FixBatteryEIS(class LAMMPS *, int, char **);
  ~FixBatteryEIS();
  void post_create();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void pre_force(int);
  void post_force(int);
  double compute_scalar();
  double compute_vector(int);

 private:
  // SOR parameters
  double omega;
  double tolerance;
  int max_iterations;
  int current_iteration;
  double convergence_error;

  // Physical constants
  double R;   // Gas constant
  double T;   // Temperature
  double F;   // Faraday constant
  double alpha_a, alpha_c;  // Transfer coefficients (not used in CC mode)

  // Material properties
  double sigma_el;  // SE ionic conductivity
  double sigma_Li;  // Li conductivity

  // Particle types
  int SE_type;
  int Li_Ctype;   // Li at cathode (bottom)
  int Li_Atype;   // Li at anode (top)
  int CC_Ctype;   // Current collector (bottom) - REFERENCE
  int CC_Atype;   // Current collector (top)

  // Boundary conditions
  double phi_ref;   // Reference potential at CC_Ctype
  double cur_app;   // Applied current density at anode (A/m²)

  // NEW: Current conservation variables
  double total_current;        // Total current through cell (A)
  double global_area_anode;    // Total contact area at anode/SE interface (m²)
  double global_area_cathode;  // Total contact area at cathode/SE interface (m²)
  double i_density_anode;      // Current density at anode interface (A/m²)
  double i_density_cathode;    // Current density at cathode interface (A/m²)

  // Per-atom property arrays
  double *phi_el;
  double *phi_el_old;
  double *phi_ed;
  double *phi_ed_old;
  double *current_Li_SE;
  double *hydrostatic_stress;
  double *contact_area_anode;    // NEW: per-particle area at anode interface
  double *contact_area_cathode;  // NEW: per-particle area at cathode interface

  // Fix property/atom handles
  class FixPropertyAtom *fix_phi_el;
  class FixPropertyAtom *fix_phi_el_old;
  class FixPropertyAtom *fix_phi_ed;
  class FixPropertyAtom *fix_phi_ed_old;
  class FixPropertyAtom *fix_current_Li_SE;
  class FixPropertyAtom *fix_hydrostatic_stress;
  class FixPropertyAtom *fix_init_flag;
  class FixPropertyAtom *fix_contact_area_anode;    // NEW
  class FixPropertyAtom *fix_contact_area_cathode;  // NEW

  // Neighbor list
  class NeighList *list;

  // State flag
  bool first_run;

  // Helper functions
  void updatePtrs();
  void solve_eis_iteration();
  void apply_reference_potential();
  void calculate_hydrostatic_stress();
  double calculate_contact_area(int i, int j);
  double calculate_current_Li_SE(int i_Li, int j_SE, double phi_ed, double phi_el, double sigma_m);
  double check_convergence();
  
  // NEW: Calculate interface currents with conservation
  void calculate_interface_currents();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix battery/eis command

Self-explanatory. Check the input script syntax and compare to the
documentation for the command.

E: Fix battery/eis requires newton pair off

This fix requires newton pair off for proper force calculation.

E: Could not find required property/atom fixes

Internal error - missing required property/atom fixes.

E: Invalid particle types for battery/eis

SE_type and Li_type must be valid atom types.

E: Invalid boundary condition types for battery/eis

Li_Ctype and Li_Atype must be valid atom types.

W: Battery EIS did not converge

The EIS solver did not converge within the maximum iterations.

W: Numerical instability detected

Numerical instability detected in potential calculation.

*/