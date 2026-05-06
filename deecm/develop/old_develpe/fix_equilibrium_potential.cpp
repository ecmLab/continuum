/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ██║   ███████║
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

    Created: Joseph Vazquez Mercado, RIT 2025
    
    Copyright 2024-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include "fix_equilibrium_potential.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "error.h"
#include "force.h"
#include "fix_property_atom.h"
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEquilibriumPotential::FixEquilibriumPotential(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  U_eq0(0.35),        // V
  R(8.31),            // J/mol·K
  T(303.0),           // K
  F(96485.0),         // C/mol
  A_param(20.0),
  B_param(-15.0),
  x_Li_max(1.00),     // Maximum Li/Si ratio
  lithium_content(NULL),
  equilibrium_potential(NULL),
  fix_lithium_content(NULL),
  fix_equilibrium_potential(NULL)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix equilibrium_potential command");

  // Parse optional arguments
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"temperature") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix equilibrium_potential command");
      T = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"U_eq0") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix equilibrium_potential command");
      U_eq0 = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"A") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix equilibrium_potential command");
      A_param = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"B") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix equilibrium_potential command");
      B_param = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"x_max") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix equilibrium_potential command");
      x_Li_max = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix equilibrium_potential command");
  }

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;
}

/* ---------------------------------------------------------------------- */

FixEquilibriumPotential::~FixEquilibriumPotential()
{
}

/* ---------------------------------------------------------------------- */

void FixEquilibriumPotential::post_create()
{
  // Register property/atom for equilibrium potential
  fix_equilibrium_potential = static_cast<FixPropertyAtom*>(modify->find_fix_property("equilibriumPotential","property/atom","scalar",0,0,style,false));
  if(!fix_equilibrium_potential) {
    const char* fixarg[10];
    fixarg[0]="equilibriumPotential";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="equilibriumPotential";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_equilibrium_potential = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }
}

/* ---------------------------------------------------------------------- */

int FixEquilibriumPotential::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEquilibriumPotential::init()
{
  // Find required fixes
  fix_lithium_content = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumContent","property/atom","scalar",0,0,style));
  if(!fix_lithium_content)
    error->all(FLERR,"Fix equilibrium_potential requires lithiumContent property");
    
  fix_equilibrium_potential = static_cast<FixPropertyAtom*>(modify->find_fix_property("equilibriumPotential","property/atom","scalar",0,0,style));
  if(!fix_equilibrium_potential)
    error->all(FLERR,"Could not find equilibriumPotential property");
    
  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixEquilibriumPotential::setup(int vflag)
{
  updatePtrs();
  // Calculate initial equilibrium potential
  calculate_equilibrium_potential();
}

/* ---------------------------------------------------------------------- */

void FixEquilibriumPotential::pre_force(int vflag)
{
  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixEquilibriumPotential::post_force(int vflag)
{
  calculate_equilibrium_potential();
}

/* ---------------------------------------------------------------------- */

void FixEquilibriumPotential::updatePtrs()
{
  lithium_content = fix_lithium_content->vector_atom;
  equilibrium_potential = fix_equilibrium_potential->vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixEquilibriumPotential::calculate_equilibrium_potential()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  
  // OLD EQUATION - COMMENTED OUT FOR REFERENCE
  // double RT_over_F = (R * T) / F;
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && type[i] == 1) {  // Only for AM particles (type 1)
      double x_li = lithium_content[i];
      
      // Prevent numerical issues
      // if (x_li <= 0.0) x_li = 1e-6;
      // if (x_li >= x_Li_max) x_li = x_Li_max - 1e-6;
      
      // Calculate normalized lithium content
      double x_norm = x_li / x_Li_max;
      
      // NEW POLYNOMIAL EQUATION For Graphite Potential vs Li Ratio
      // U_eq = 2.60481*10^4*x_norm^12 - 1.65849*10^5*x_norm^11 + 4.6562*10^5*x_norm^10 
      //        - 7.59309*10^5*x_norm^9 + 7.97656*10^5*x_norm^8 - 5.6595*10^5*x_norm^7 
      //        + 2.76551*10^5*x_norm^6 - 9.31284*10^4*x_norm^5 + 2.12529*10^4*x_norm^4 
      //        - 3.16719*10^3*x_norm^3 + 2.89822*10^2*x_norm^2 - 1.5121*10^1*x_norm 
      //        + 5.5141*10^-1
      
      double x_norm2 = x_norm * x_norm;
      double x_norm3 = x_norm2 * x_norm;
      double x_norm4 = x_norm3 * x_norm;
      // double x_norm5 = x_norm4 * x_norm;
      // double x_norm6 = x_norm5 * x_norm;
      // double x_norm7 = x_norm6 * x_norm;
      // double x_norm8 = x_norm7 * x_norm;
      // double x_norm9 = x_norm8 * x_norm;
      // double x_norm10 = x_norm9 * x_norm;
      // double x_norm11 = x_norm10 * x_norm;
      // double x_norm12 = x_norm11 * x_norm;
      
      double numerator1 = -17709.3425099138 * x_norm4 + 
                       49563.5287557803 * x_norm3 + 
                       -84640.8295904228 * x_norm2 + 
                       -4980.8894321468 * x_norm + 
                       56960.9582619908;

      double denominator1 = x_norm3 + 
                         -14391.9820106230 * x_norm2 + 
                         1641.3022786177 * x_norm + 
                         12525.8318098236;
      
      equilibrium_potential[i] = numerator1 / denominator1;

      // equilibrium_potential[i] = 2.60481e4 * x_norm12 - 1.65849e5 * x_norm11 + 4.6562e5 * x_norm10
      //                           - 7.59309e5 * x_norm9 + 7.97656e5 * x_norm8 - 5.6595e5 * x_norm7
      //                           + 2.76551e5 * x_norm6 - 9.31284e4 * x_norm5 + 2.12529e4 * x_norm4
      //                           - 3.16719e3 * x_norm3 + 2.89822e2 * x_norm2 - 1.5121e1 * x_norm
      //                           + 5.5141e-1;
      
      // OLD EQUATION - COMMENTED OUT FOR REFERENCE
      // double term1 = log((1.0 - x_li/x_Li_max) / (x_li/x_Li_max));
      // double numerator = (B_param - 1.0) * x_li * x_li / (x_Li_max * x_Li_max) + 2.0 * x_li / x_Li_max - 1.0;
      // double denominator = 1.0 + (B_param - 1.0) * x_li / x_Li_max;
      // denominator = denominator * denominator;
      //
      // if (fabs(denominator) < 1e-10) {
      //   error->warning(FLERR,"Near-zero denominator in equilibrium potential calculation");
      //   denominator = 1e-10;
      // }
      //
      // double term2 = A_param * numerator / denominator;
      // 
      // equilibrium_potential[i] = U_eq0 + RT_over_F * (term1 + term2);
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixEquilibriumPotential::compute_scalar()
{
  // Return average equilibrium potential for AM particles
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double sum = 0.0;
  int count = 0;
  
  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit && type[i] == 1) {
      sum += equilibrium_potential[i];
      count++;
    }
  }
  
  double all_sum, all_count;
  MPI_Allreduce(&sum,&all_sum,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&count,&all_count,1,MPI_INT,MPI_SUM,world);
  
  if (all_count > 0) return all_sum/all_count;
  return 0.0;
}
