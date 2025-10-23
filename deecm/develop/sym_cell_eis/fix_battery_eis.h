/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    
    Copyright 2024-     DCS Computing GmbH, Linz
    
    Modified for NMC cathode simulation with dual potential solving
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
  virtual ~FixBatteryEIS();
  virtual void post_create();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void pre_force(int);
  void post_force(int);
  double compute_scalar();
  double compute_vector(int);

  virtual void updatePtrs();
  
 protected:
  // EIS parameters
  double omega;              // EIS relaxation parameter (1.9)
  double tolerance;          // Convergence tolerance
  int max_iterations;        // Maximum EIS iterations
  int current_iteration;     // Current iteration count
  double convergence_error;  // Current convergence error
  
  // Physical parameters
  double R;                  // Gas constant (8.31 J/mol·K)
  double T;                  // Temperature (303 K)
  double F;                  // Faraday constant (96485 C/mol)
  
  // Conductivity parameters
  double sigma_el;           // Electrolyte conductivity (SE) (S/m)
  double sigma_ed_SE;        // Electronic conductivity for CBD/SE (S/m)
  double sigma_ed_CC;        // Electronic conductivity for CC (S/m)
  
  double alpha_a;            // Anodic transfer coefficient (0.5)
  double alpha_c;            // Cathodic transfer coefficient (0.5)
  
  // Boundary conditions
  double phi_el_BC_bottom;   // Electrolyte potential at bottom boundary
  double phi_el_BC_top;      // Electrolyte potential at top boundary
  double phi_ed_BC_anode;    // Electronic potential at anode (0V)
  
  // Property pointers - Electrolyte
  double *phi_el;            // Electrolyte potential
  double *phi_el_old;        // Old electrolyte potential (for convergence)
  
  // Property pointers - Electronic
  double *phi_ed;            // Electronic potential
  double *phi_ed_old;        // Old electronic potential (for convergence)
  
  // Property pointers - Other
  double *current_Li_SE;     // Current density from Li to SE
  double *hydrostatic_stress; // Hydrostatic stress on Li particles
  
  // Fix pointers - Electrolyte
  class FixPropertyAtom *fix_phi_el;
  class FixPropertyAtom *fix_phi_el_old;
  
  // Fix pointers - Electronic
  class FixPropertyAtom *fix_phi_ed;
  class FixPropertyAtom *fix_phi_ed_old;
  
  // Fix pointers - Other
  class FixPropertyAtom *fix_current_Li_SE;
  class FixPropertyAtom *fix_hydrostatic_stress;
  class FixPropertyAtom *fix_init_flag;
  
  // Particle type groups
  int BC_bottom_type;        // Atom type for bottom BC particles (CC)
  int BC_top_type;           // Atom type for top BC particles (Anode)
  int groupbit_SE;           // Group bit for SE particles
  int SE_type;               // Atom type for SE particles
  
  // Neighbor list
  class NeighList *list;
  
  // Methods
  void solve_eis_iteration();
  void apply_boundary_conditions();
  void calculate_hydrostatic_stress();
  double calculate_contact_area(int, int);
  double calculate_current_Li_SE(int, int, double, double, double);
  double check_convergence();
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

BC_bottom_type and BC_top_type must be valid atom types.

W: Battery EIS did not converge

The EIS solver did not converge within the maximum iterations.

W: Numerical instability detected

Numerical instability detected in potential calculation.

*/