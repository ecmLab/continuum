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
  // Physical constants
  double F;                  // Faraday constant (96485 C/mol)

  // User-configurable material parameters
  double c_li_max;           // Maximum Li concentration (mol/m³)
  double D_Li;               // Li self-diffusion coefficient (m²/s)
  double Omega_Li;           // Molar volume of lithium (m³/mol)

  // User-configurable simulation parameters
  double s_factor;           // Electrochemical time scaling factor
  double time_conv;          // Time unit conversion factor to seconds

  // Content clamping (optional)
  bool clamp_content;        // Whether to clamp lithium content
  double clamp_min;          // Minimum lithium content bound
  double clamp_max;          // Maximum lithium content bound

  // Lithium content parameters (retrieved from fix_property_atom_lithium_content)
  double initial_lithium_content;
  double target_lithium_content;
  double max_lithium_content;

  // Property pointers
  double *lithium_content;
  double *lithium_concentration;
  double *current_Li_SE;
  double *diffusion_coefficient;
  double *lithium_flux;
  double *li_mols;
  
  // Fix pointers
  class FixPropertyAtom *fix_lithium_content;
  class FixPropertyAtom *fix_lithium_concentration;
  class FixPropertyAtom *fix_current_Li_SE;
  class FixPropertyAtom *fix_diffusion_coefficient;
  class FixPropertyAtom *fix_lithium_flux;
  class FixPropertyAtom *fix_li_mols;
  class FixPropertyAtomLithiumContent *fix_lithium_content_manager;
  
  // Particle type
  int Li_Atype;
  int Li_Ctype;
  
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

E: Illegal fix lithium_diffusion command: Li_type requires 2 values

The Li_type keyword requires anode_type and cathode_type.

E: Illegal fix lithium_diffusion command: unknown keyword

An unrecognized keyword was used.

E: D_Li must be positive

The diffusion coefficient must be a positive number.

E: Omega_Li must be positive

The molar volume must be a positive number.

E: time_conv must be positive

The time conversion factor must be a positive number.

E: clamp_content min must be less than max

The minimum clamp value must be strictly less than the maximum.

E: Fix lithium_diffusion requires lithiumContent property

Fix property/atom/lithium_content must be defined before fix lithium_diffusion.

E: Fix lithium_diffusion requires fix property/atom/lithium_content

Fix property/atom/lithium_content must be defined before fix lithium_diffusion.

*/