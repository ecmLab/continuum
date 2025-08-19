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
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(battery/sor,FixBatterySOR)

#else

#ifndef LMP_FIX_BATTERY_SOR_H
#define LMP_FIX_BATTERY_SOR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBatterySOR : public Fix {
 public:
  FixBatterySOR(class LAMMPS *, int, char **);
  virtual ~FixBatterySOR();
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
  // SOR parameters
  double omega;              // SOR relaxation parameter (1.9)
  double tolerance;          // Convergence tolerance
  int max_iterations;        // Maximum SOR iterations
  int current_iteration;     // Current iteration count
  double convergence_error;  // Current convergence error
  
  // Physical parameters
  double R;                  // Gas constant (8.31 J/mol·K)
  double T;                  // Temperature (303 K)
  double F;                  // Faraday constant (96485 C/mol)
  double sigma_SE;           // SE conductivity (0.05 S/m)
  double alpha_a;            // Anodic transfer coefficient (0.5)
  double alpha_c;            // Cathodic transfer coefficient (0.5)
  
  // Boundary conditions
  double phi_el_BC_bottom;   // Electric potential at bottom boundary (0.0 V)
  double phi_el_BC_top;      // Electric potential at top boundary (0.04 V)
  
  // Property pointers
  double *phi_el;            // Electric potential
  double *phi_el_old;        // Old electric potential (for convergence check)
  double *equilibrium_potential;
  double *exchange_current_density;
  double *current_AM_SE;     // Current density from AM to SE
  double *hydrostatic_stress; // Hydrostatic stress on AM particles
  
  // Fix pointers
  class FixPropertyAtom *fix_phi_el;
  class FixPropertyAtom *fix_phi_el_old;
  class FixPropertyAtom *fix_equilibrium_potential;
  class FixPropertyAtom *fix_exchange_current_density;
  class FixPropertyAtom *fix_current_AM_SE;
  class FixPropertyAtom *fix_hydrostatic_stress;
  
  // Particle type groups
  int BC_bottom_type;	     // Group bit for CC1 particles
  int BC_top_type;           // Group bit for CC2 particles
  int groupbit_SE;           // Group bit for SE particles
  int groupbit_AM;           // Group bit for AM particles
  int SE_type;               // Atom type for SE particles
  int AM_type;               // Atom type for AM particles
  
  // Neighbor list
  class NeighList *list;
  
  // Methods
  void solve_sor_iteration();
  void apply_boundary_conditions();
  void calculate_hydrostatic_stress();
  double calculate_contact_area(int, int);
  double calculate_current_AM_SE(int, int, double, double, double);
  double check_convergence();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix battery/sor command

Self-explanatory. Check the input script syntax and compare to the
documentation for the command.

E: Fix battery/sor requires newton pair off

This fix requires newton pair off for proper force calculation.

E: Fix battery/sor requires fix equilibrium_potential

Fix equilibrium_potential must be defined before fix battery/sor.

E: Fix battery/sor requires fix exchange_current_density

Fix exchange_current_density must be defined before fix battery/sor.

E: Could not find required property/atom fixes

Internal error - missing required property/atom fixes.

E: Invalid particle types for battery/sor

SE_type and AM_type must be valid atom types.

W: Battery SOR did not converge

The SOR solver did not converge within the maximum iterations.

*/
