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

FixStyle(equilibrium_potential,FixEquilibriumPotential)

#else

#ifndef LMP_FIX_EQUILIBRIUM_POTENTIAL_H
#define LMP_FIX_EQUILIBRIUM_POTENTIAL_H

#include "fix.h"

namespace LAMMPS_NS {

// Material type enumeration
enum MaterialType {
  MATERIAL_AG = 0,
  MATERIAL_C = 1,
  MATERIAL_NMC811 = 2
};

class FixEquilibriumPotential : public Fix {
 public:
  FixEquilibriumPotential(class LAMMPS *, int, char **);
  virtual ~FixEquilibriumPotential();
  virtual void post_create();
  int setmask();
  void init();
  void setup(int);
  void pre_force(int vflag);
  void post_force(int vflag);
  double compute_scalar();

  virtual void updatePtrs();
  
 protected:
  // Model parameters from the flowchart
  double R;          // Gas constant (8.31 J/mol·K)
  double T;          // Temperature (303 K)
  double F;          // Faraday constant (96485 C/mol)
  double x_Li_max;   // Maximum Li/Si ratio (3.75)
  
  // Material type selection
  int material_type;
  
  // Property pointers
  double *lithium_content;
  double *equilibrium_potential;
  
  // Fix pointers
  class FixPropertyAtom *fix_lithium_content;
  class FixPropertyAtom *fix_equilibrium_potential;
  
  void calculate_equilibrium_potential();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix equilibrium_potential command

Self-explanatory. Check the input script syntax and compare to the
documentation for the command.

E: Fix equilibrium_potential requires fix property/atom/lithium_content

Fix property/atom/lithium_content must be defined before fix equilibrium_potential.

E: Could not find fix property/atom equilibriumPotential

Internal error.

E: Could not find fix property/atom lithiumContent

Internal error.

E: Unknown material type for fix equilibrium_potential

Valid material types are: Ag, C, NMC811

*/