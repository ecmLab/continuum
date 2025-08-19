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
  // Diffusion parameters from equations
  double D_min;              // Minimum diffusion coefficient (1e-15 m²/s)
  double D_poor;             // Poor diffusion coefficient (1e-14 m²/s)
  double D_rich;             // Rich diffusion coefficient (1e-13 m²/s)
  double kappa;              // Exponential parameter (-1)
  double F;                  // Faraday constant (96485 C/mol)
  double c_li_max;           // Maximum Li concentration (83874 mol/m³)

  // Lithium content parameters (retrieved from fix_property_atom_lithium_content)
  double initial_lithium_content;  // Initial Li/Si ratio
  double target_lithium_content;   // Target Li/Si ratio for charging
  double max_lithium_content;      // Maximum Li/Si ratio

  // Property pointers
  double *lithium_content;
  double *lithium_concentration;
  double *current_AM_SE;
  double *diffusion_coefficient;
  double *lithium_flux;
  double *silicon_content;
  
  // Fix pointers
  class FixPropertyAtom *fix_lithium_content;
  class FixPropertyAtom *fix_lithium_concentration;
  class FixPropertyAtom *fix_current_AM_SE;
  class FixPropertyAtom *fix_diffusion_coefficient;
  class FixPropertyAtom *fix_lithium_flux;
  class FixPropertyAtom *fix_silicon_content;
  class FixPropertyAtomLithiumContent *fix_lithium_content_manager;
  
  // Particle type
  int AM_type;
  
  // Neighbor list
  class NeighList *list;
  
  // Methods
  void calculate_diffusion_coefficient();
  void update_lithium_content();
  double calculate_contact_area(int, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix lithium_diffusion command

Self-explanatory. Check the input script syntax and compare to the
documentation for the command.

E: Fix lithium_diffusion requires fix property/atom/lithium_content

Fix property/atom/lithium_content must be defined before fix lithium_diffusion.

E: Fix lithium_diffusion requires fix exchange_current_density

Fix exchange_current_density must be defined before fix lithium_diffusion.

E: Fix lithium_diffusion requires fix battery/sor

Fix battery/sor must be defined before fix lithium_diffusion.

*/
