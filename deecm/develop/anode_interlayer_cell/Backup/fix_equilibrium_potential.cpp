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
  R(8.31),            // J/mol·K
  T(303.0),           // K
  F(96485.0),         // C/mol
  x_Li_max(1.00),     // Maximum Li/Si ratio
  material_type(MATERIAL_AG),  // Default to Ag
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
    } else if (strcmp(arg[iarg],"x_max") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix equilibrium_potential command");
      x_Li_max = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"material") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix equilibrium_potential command");
      if (strcmp(arg[iarg+1],"Ag") == 0) {
        material_type = MATERIAL_AG;
      } else if (strcmp(arg[iarg+1],"C") == 0) {
        material_type = MATERIAL_C;
      } else if (strcmp(arg[iarg+1],"NMC811") == 0) {
        material_type = MATERIAL_NMC811;
      } else {
        error->all(FLERR,"Unknown material type for fix equilibrium_potential. Valid types: Ag, C, NMC811");
      }
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
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && type[i] == 1) {  // Only for AM particles (type 1)
      double x_li = lithium_content[i];
      double x_norm = x_li / x_Li_max;
      
      double x_norm2 = x_norm * x_norm;
      double x_norm3 = x_norm2 * x_norm;
      double x_norm4 = x_norm3 * x_norm;
      
      double numerator, denominator;
      
      if (material_type == MATERIAL_AG) {
        // Ag: V = (-0.1631924300*x^4 + 0.4016544523*x^3 + -0.3062847559*x^2 + 0.0748140933*x + 0.0000767930)/(x^3 + -1.1531920913*x^2 + 0.3387447153*x + -0.0008769410)
        numerator = -0.1631924300 * x_norm4 + 
                     0.4016544523 * x_norm3 + 
                    -0.3062847559 * x_norm2 + 
                     0.0748140933 * x_norm + 
                     0.0000767930;
        denominator = x_norm3 + 
                     -1.1531920913 * x_norm2 + 
                      0.3387447153 * x_norm + 
                     -0.0008769410;
      } else if (material_type == MATERIAL_C) {
        // C: V = (-574.4685917325*x^4 + 438.1745026465*x^3 + -122.1174477257*x^2 + 271.4976114422*x + 10.8530532232)/(x^3 + 262.5473706486*x^2 + 2710.3782301839*x + -117.7939698484)
        numerator = -574.4685917325 * x_norm4 + 
                     438.1745026465 * x_norm3 + 
                    -122.1174477257 * x_norm2 + 
                     271.4976114422 * x_norm + 
                      10.8530532232;
        denominator = x_norm3 + 
                      262.5473706486 * x_norm2 + 
                     2710.3782301839 * x_norm + 
                     -117.7939698484;
      } else {  // MATERIAL_NMC811
        // NMC811: V = (-17709.3425099138*x^4 + 49563.5287557803*x^3 + -84640.8295904228*x^2 + -4980.8894321468*x + 56960.9582619908)/(x^3 + -14391.9820106230*x^2 + 1641.3022786177*x + 12525.8318098236)
        numerator = -17709.3425099138 * x_norm4 + 
                     49563.5287557803 * x_norm3 + 
                    -84640.8295904228 * x_norm2 + 
                     -4980.8894321468 * x_norm + 
                     56960.9582619908;
        denominator = x_norm3 + 
                     -14391.9820106230 * x_norm2 + 
                       1641.3022786177 * x_norm + 
                      12525.8318098236;
      }
      
      equilibrium_potential[i] = numerator / denominator;
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