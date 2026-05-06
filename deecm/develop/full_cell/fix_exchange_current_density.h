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

FixStyle(exchange_current_density,FixExchangeCurrentDensity)

#else

#ifndef LMP_FIX_EXCHANGE_CURRENT_DENSITY_H
#define LMP_FIX_EXCHANGE_CURRENT_DENSITY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixExchangeCurrentDensity : public Fix {
 public:
  FixExchangeCurrentDensity(class LAMMPS *, int, char **);
  virtual ~FixExchangeCurrentDensity();
  virtual void post_create();
  int setmask();
  void init();
  void setup(int);              // ADD THIS LINE
  void pre_force(int vflag);
  void post_force(int vflag);
  double compute_scalar();

  virtual void updatePtrs();
  
 protected:
  // Model parameters from the flowchart
  double F;              // Faraday constant (96485 C/mol)
  double k_r;            // Rate constant (6e-12 m/s)
  double c_li_max;       // Maximum Li concentration (83874 mol/m³)
  double c_Li_plus_q;    // Li+ concentration in electrolyte (50.608 mol/m³)
  double alpha_a;        // Anodic transfer coefficient (0.5)
  double alpha_c;        // Cathodic transfer coefficient (0.5)
  double x_Li_max;       // Maximum Li/Si ratio (3.75)
  
  // Property pointers
  double *lithium_content;
  double *exchange_current_density;
  double *lithium_concentration;
  
  // Fix pointers
  class FixPropertyAtom *fix_lithium_content;
  class FixPropertyAtom *fix_exchange_current_density;
  class FixPropertyAtom *fix_lithium_concentration;
  
  void calculate_exchange_current_density();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix exchange_current_density command

Self-explanatory. Check the input script syntax and compare to the
documentation for the command.

E: Fix exchange_current_density requires fix property/atom/lithium_content

Fix property/atom/lithium_content must be defined before fix exchange_current_density.

E: Could not find fix property/atom exchangeCurrentDensity

Internal error.

E: Could not find fix property/atom lithiumContent

Internal error.

E: Could not find fix property/atom lithiumConcentration

Internal error.

*/
