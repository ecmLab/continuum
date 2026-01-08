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

FixStyle(potential/sor,FixPotentialSOR)

#else

#ifndef LMP_FIX_POTENTIAL_SOR_H
#define LMP_FIX_POTENTIAL_SOR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPotentialSOR : public Fix {
 public:
  FixPotentialSOR(class LAMMPS *, int, char **);
  virtual ~FixPotentialSOR();
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
  int max_iterations;        // Maximum EIS iterations
  int current_iteration;     // Current iteration count
  double convergence_error;  // Current convergence error
  
  // Physical parameters
  double R;                  // Gas constant (8.31 J/mol·K)
  double T;                  // Temperature (303 K)
  double F;                  // Faraday constant (96485 C/mol)
  
  // Conductivity parameters
  double sigma_el;           // Electrolyte conductivity (SE) (S/m)
  double sigma_Li;  // Li conductivity
  double sigma_ed_SE;        // Electronic conductivity for CBD/SE (S/m)
  double sigma_ed_CC;        // Electronic conductivity for CC (S/m)
  
  double alpha_a;            // Anodic transfer coefficient (0.5)
  double alpha_c;            // Cathodic transfer coefficient (0.5)
  
  // Boundary conditions
  double phi_el_BC_Cat;   // Electrolyte potential at bottom boundary
  double phi_el_BC_An;      // Electrolyte potential at top boundary
  double phi_ed_BC_anode;    // Electronic potential at anode (0V)
  double cur_app;          // Applied current density (A/m2)
  
  // Property pointers - Electrolyte
  double *phi_el;            // Electrolyte potential
  double *phi_el_old;        // Old electrolyte potential (for convergence)
  
  // Property pointers - Electronic
  double *phi_ed;            // Electronic potential
  double *phi_ed_old;        // Old electronic potential (for convergence)
  
  // Property pointers - Other
  double *current_Li_SE;     // Current density from Li to SE
  double *hydrostatic_stress; // Hydrostatic stress on Li particles

  // Current distribution variables
  double total_current;
  double global_area_anode;
  double global_area_cathode;
  double i_density_anode;
  double i_density_cathode;

  // Per-particle contact area pointers
  double *contact_area_anode;
  double *contact_area_cathode;
  
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
  
  // Fix pointers for contact areas
  class FixPropertyAtom *fix_contact_area_anode;
  class FixPropertyAtom *fix_contact_area_cathode;

  // Particle type groups
  int Li_Ctype;        // Atom type for bottom BC particles (CC)
  int Li_Atype;           // Atom type for top BC particles (Anode)
  int groupbit_SE;           // Group bit for SE particles
  int SE_type;               // Atom type for SE particles
  
  // State flag
  bool first_run;

  // Neighbor list
  class NeighList *list;
  
  // Methods
  void solve_eis_iteration();
  void apply_boundary_conditions();
  void calculate_hydrostatic_stress();
  double calculate_contact_area(int, int);
  double calculate_current_Li_SE(int, int, double, double, double);
  double check_convergence();
  // NEW: Calculate interface currents with conservation
  void calculate_interface_currents();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix potential/sor command

Self-explanatory. Check the input script syntax and compare to the
documentation for the command.

E: Fix potential/sor requires newton pair off

This fix requires newton pair off for proper force calculation.

E: Could not find required property/atom fixes

Internal error - missing required property/atom fixes.

E: Invalid particle types for potential/sor

SE_type and Li_type must be valid atom types.

E: Invalid boundary condition types for potential/sor

Li_Ctype and Li_Atype must be valid atom types.

W: Battery EIS did not converge

The EIS solver did not converge within the maximum iterations.

W: Numerical instability detected

Numerical instability detected in potential calculation.

*/