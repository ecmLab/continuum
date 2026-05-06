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

#include "fix_exchange_current_density.h"
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

FixExchangeCurrentDensity::FixExchangeCurrentDensity(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  F(96485.0),              // C/mol
  k_r(6e-12),              // m/s
  c_li_max(83874.0),       // mol/m³
  c_Li_plus_q(50608.0),     // mol/m³
  alpha_a(0.5),
  alpha_c(0.5),
  x_Li_max(1.00),         // Maximum Li/Si ratio
  lithium_content(NULL),
  exchange_current_density(NULL),
  lithium_concentration(NULL),
  fix_lithium_content(NULL),
  fix_exchange_current_density(NULL),
  fix_lithium_concentration(NULL)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix exchange_current_density command");

  // Parse optional arguments
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"k_r") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix exchange_current_density command");
      k_r = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"c_li_max") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix exchange_current_density command");
      c_li_max = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"c_electrolyte") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix exchange_current_density command");
      c_Li_plus_q = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"alpha_a") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix exchange_current_density command");
      alpha_a = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"alpha_c") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix exchange_current_density command");
      alpha_c = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"x_max") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix exchange_current_density command");
      x_Li_max = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix exchange_current_density command");
  }

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;
}

/* ---------------------------------------------------------------------- */

FixExchangeCurrentDensity::~FixExchangeCurrentDensity()
{
}

/* ---------------------------------------------------------------------- */

void FixExchangeCurrentDensity::post_create()
{
  // Register property/atom for exchange current density
  fix_exchange_current_density = static_cast<FixPropertyAtom*>(modify->find_fix_property("exchangeCurrentDensity","property/atom","scalar",0,0,style,false));
  if(!fix_exchange_current_density) {
    const char* fixarg[10];
    fixarg[0]="exchangeCurrentDensity";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="exchangeCurrentDensity";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_exchange_current_density = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  // Register property/atom for lithium concentration
  fix_lithium_concentration = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumConcentration","property/atom","scalar",0,0,style,false));
  if(!fix_lithium_concentration) {
    const char* fixarg[10];
    fixarg[0]="lithiumConcentration";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="lithiumConcentration";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_lithium_concentration = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }
}

/* ---------------------------------------------------------------------- */

int FixExchangeCurrentDensity::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixExchangeCurrentDensity::init()
{
  // Find required fixes
  fix_lithium_content = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumContent","property/atom","scalar",0,0,style));
  if(!fix_lithium_content)
    error->all(FLERR,"Fix exchange_current_density requires lithiumContent property");
    
  fix_exchange_current_density = static_cast<FixPropertyAtom*>(modify->find_fix_property("exchangeCurrentDensity","property/atom","scalar",0,0,style));
  if(!fix_exchange_current_density)
    error->all(FLERR,"Could not find exchangeCurrentDensity property");
    
  fix_lithium_concentration = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumConcentration","property/atom","scalar",0,0,style));
  if(!fix_lithium_concentration)
    error->all(FLERR,"Could not find lithiumConcentration property");
    
  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixExchangeCurrentDensity::setup(int vflag)
{
  updatePtrs();
  // Calculate initial exchange current density
  calculate_exchange_current_density();
}

/* ---------------------------------------------------------------------- */

void FixExchangeCurrentDensity::pre_force(int vflag)
{
  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixExchangeCurrentDensity::post_force(int vflag)
{
  calculate_exchange_current_density();
}

/* ---------------------------------------------------------------------- */

void FixExchangeCurrentDensity::updatePtrs()
{
  lithium_content = fix_lithium_content->vector_atom;
  exchange_current_density = fix_exchange_current_density->vector_atom;
  lithium_concentration = fix_lithium_concentration->vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixExchangeCurrentDensity::calculate_exchange_current_density()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && type[i] == 1) {  // Only for AM particles (type 1)
      double x_li = lithium_content[i];
      
      // Calculate lithium concentration c_Li,p from molar ratio
      lithium_concentration[i] = x_li * c_li_max / x_Li_max;
      double c_li = lithium_concentration[i];
      
      // Prevent numerical issues
      // if (c_li <= 0.0) c_li = 1e-6;
      // if (c_li >= c_li_max) c_li = c_li_max - 1e-6;
      
      // Calculate exchange current density using the formula from the flowchart
      // i_p,q^(0) = F * k_r * (c_Li,p / c_Li,p^(max) * c_Li+,q)^alpha_a * (c_Li,p)^(1-alpha_c)
      // double term1 = pow((c_li / c_li_max) * c_Li_plus_q, alpha_a);
      // double term2 = pow(c_li, 1.0 - alpha_c);
      
      // exchange_current_density[i] = F * k_r * term1 * term2; // A/m²
      exchange_current_density[i] = 260; // A/m² - Ag Constant Value (Nuhu Tafel Plot Exp. 11/17/2025, Other ref is https://doi.org/10.1002/adma.202303489)
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixExchangeCurrentDensity::compute_scalar()
{
  // Return average exchange current density for AM particles
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double sum = 0.0;
  int count = 0;
  
  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit && type[i] == 1) {
      sum += exchange_current_density[i];
      count++;
    }
  }
  
  double all_sum, all_count;
  MPI_Allreduce(&sum,&all_sum,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&count,&all_count,1,MPI_INT,MPI_SUM,world);
  
  if (all_count > 0) return all_sum/all_count;
  return 0.0;
}
