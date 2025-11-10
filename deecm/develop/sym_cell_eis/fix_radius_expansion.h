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
  // Volume parameters from equations
  double Omega_Li;           // Molar volume of Li (12.97e-6 m³/mol)
  double Omega_Li_eff;       // Effective molar volume of Li in Li (12.97e-6 m³/mol)
  
  // Update control
  int update_frequency;      // Steps between radius updates
  int max_updates;           // Maximum number of radius updates (10)
  int update_count;          // Current update count
  int next_update_step;      // Next step to update radius
  bool manual_trigger;       // Manual trigger for radius update
  
  // Property pointers
  double *lithium_content;
  double *li_mols;      // Moles of lithium (from lithium_diffusion)
  double *initial_radius;
  double *particle_volume;
  
  // Fix pointers
  class FixPropertyAtom *fix_lithium_content;
  class FixPropertyAtom *fix_li_mols;
  class FixPropertyAtom *fix_initial_radius;
  class FixPropertyAtom *fix_particle_volume;
  
  // Particle type
  int Li_Atype;
  int Li_Ctype;
  
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

E: Fix radius_expansion requires fix property/atom/lithium_content

Fix property/atom/lithium_content must be defined before fix radius_expansion.

E: Maximum radius updates reached

The fix has performed the maximum allowed number of radius updates.

*/
