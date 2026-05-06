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

FixStyle(property/atom/lithium_content,FixPropertyAtomLithiumContent)

#else

#ifndef LMP_FIX_PROPERTY_ATOM_LITHIUM_CONTENT_H
#define LMP_FIX_PROPERTY_ATOM_LITHIUM_CONTENT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPropertyAtomLithiumContent : public Fix {
 public:
  FixPropertyAtomLithiumContent(class LAMMPS *, int, char **);
  virtual ~FixPropertyAtomLithiumContent();
  virtual void post_create();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  double compute_scalar();
  double compute_vector(int);

  virtual void updatePtrs();
  double* get_lithium_content() { return lithium_content; }
  double get_max_lithium_content() { return max_lithium_content; }
  
 private:
  double initial_lithium_content;   // Starting Li/Si ratio (0.075)
  double target_lithium_content;    // End Li/Si ratio (3.3)
  double max_lithium_content;       // Maximum Li/Si ratio (3.75)
  
  double *lithium_content;          // Per-atom Li/Si molar ratio
  class FixPropertyAtom *fix_lithium_content;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix property/atom/lithium_content command

Self-explanatory. Check the input script syntax and compare to the
documentation for the command.

E: Invalid initial lithium content

Initial lithium content must be between 0 and maximum lithium content.

E: Invalid target lithium content

Target lithium content must be between 0 and maximum lithium content.

E: Invalid maximum lithium content

Maximum lithium content must be greater than 0.

E: Could not find fix property/atom lithiumContent

Internal error.

*/
