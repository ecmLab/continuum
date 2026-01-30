/* ----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine
    Contributing author: Joseph Vazquez Mercado, RIT 2025
    Copyright 2024-     DCS Computing GmbH, Linz
    
    Updated: Support for multiple particle types (AM and CB) with 
    different OCV formulas
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
  R(8.314),           // J/(mol·K)
  T(303.0),           // K
  F(96485.0),         // C/mol
  x_Li_max_AM(1.0),   // Max Li content for AM
  x_Li_max_CB(1.0),   // Max Li content for CB
  U_OCV_LM(0.0),      // OCV for lithium metal (0V vs Li/Li+)
  AM_material(MATERIAL_AG),
  CB_material(MATERIAL_C),
  AM_type(1),
  CB_type(2),
  LM_type(4),
  lithium_content(NULL),
  equilibrium_potential(NULL),
  fix_lithium_content(NULL),
  fix_equilibrium_potential(NULL)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix equilibrium_potential command");

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"temperature") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix equilibrium_potential command");
      T = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"x_max_AM") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix equilibrium_potential command");
      x_Li_max_AM = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"x_max_CB") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix equilibrium_potential command");
      x_Li_max_CB = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"U_OCV_LM") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix equilibrium_potential command");
      U_OCV_LM = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"AM_material") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix equilibrium_potential command");
      if (strcmp(arg[iarg+1],"Ag") == 0) {
        AM_material = MATERIAL_AG;
      } else if (strcmp(arg[iarg+1],"C") == 0) {
        AM_material = MATERIAL_C;
      } else if (strcmp(arg[iarg+1],"NMC811") == 0) {
        AM_material = MATERIAL_NMC811;
      } else {
        error->all(FLERR,"Unknown AM material type. Valid: Ag, C, NMC811");
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"CB_material") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix equilibrium_potential command");
      if (strcmp(arg[iarg+1],"Ag") == 0) {
        CB_material = MATERIAL_AG;
      } else if (strcmp(arg[iarg+1],"C") == 0) {
        CB_material = MATERIAL_C;
      } else if (strcmp(arg[iarg+1],"NMC811") == 0) {
        CB_material = MATERIAL_NMC811;
      } else {
        error->all(FLERR,"Unknown CB material type. Valid: Ag, C, NMC811");
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"AM_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix equilibrium_potential command");
      AM_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"CB_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix equilibrium_potential command");
      CB_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"LM_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix equilibrium_potential command");
      LM_type = force->inumeric(FLERR,arg[iarg+1]);
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
  fix_equilibrium_potential = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("equilibriumPotential","property/atom","scalar",0,0,style,false));
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
  fix_lithium_content = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("lithiumContent","property/atom","scalar",0,0,style));
  if(!fix_lithium_content)
    error->all(FLERR,"Fix equilibrium_potential requires lithiumContent property");
    
  fix_equilibrium_potential = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("equilibriumPotential","property/atom","scalar",0,0,style));
  if(!fix_equilibrium_potential)
    error->all(FLERR,"Could not find equilibriumPotential property");
    
  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixEquilibriumPotential::setup(int vflag)
{
  updatePtrs();
  calculate_equilibrium_potential();
  fix_equilibrium_potential->do_forward_comm();
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
  fix_equilibrium_potential->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

void FixEquilibriumPotential::updatePtrs()
{
  lithium_content = fix_lithium_content->vector_atom;
  equilibrium_potential = fix_equilibrium_potential->vector_atom;
}

/* ---------------------------------------------------------------------- */

double FixEquilibriumPotential::compute_OCV(double x_li, double x_max, int material)
{
  double x_norm = x_li / x_max;
  
  // // Clamp x_norm to valid range
  // if (x_norm < 0.001) x_norm = 0.001;
  // if (x_norm > 0.999) x_norm = 0.999;
  
  double x2 = x_norm * x_norm;
  double x3 = x2 * x_norm;
  double x4 = x3 * x_norm;
  double x5 = x4 * x_norm;
  double x6 = x5 * x_norm;
  double x7 = x6 * x_norm;
  double x8 = x7 * x_norm;
  double x9 = x8 * x_norm;
  
  double num, den;
  
  if (material == MATERIAL_AG) {
    num = 0.4980318604*x6 + -0.6468695682*x5 + 0.3082756996*x4 + -0.0595231213*x3 + 0.0051736954*x2 + -0.0001934195*x_norm + 0.0000026074;
    den = x8 + 26.7713969498*x7 + -23.2071887899*x6 + 8.2525713368*x5 + -1.1034187782*x4 + 0.0355382579*x3 + 0.0032433057*x2 + -0.0002229515*x_norm + 0.0000036652;
  } else if (material == MATERIAL_C) {
    num = -17.6001525048*x9 + 66.6428247669*x8 + -105.9285748917*x7 + 91.1122366366*x6 + -45.5078952362*x5 + 13.3337219668*x4 + -2.2539518017*x3 + 0.2176770426*x2 + -0.0111586893*x_norm + 0.0002353456;
    den = x8 + -28.1174547327*x7 + 62.3756782904*x6 + -51.5130168878*x5 + 19.8169524758*x4 + -3.7413389704*x3 + 0.3746346388*x2 + -0.0193247118*x_norm + 0.0004061816;
  } else { // MATERIAL_NMC811
    num = -17709.3425099138*x4 + 49563.5287557803*x3 - 84640.8295904228*x2 - 4980.8894321468*x_norm + 56960.9582619908;
    den = x3 - 14391.9820106230*x2 + 1641.3022786177*x_norm + 12525.8318098236;
  }
  
  return num / den;
}

/* ---------------------------------------------------------------------- */

void FixEquilibriumPotential::calculate_equilibrium_potential()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    
    if (type[i] == AM_type) {
      equilibrium_potential[i] = compute_OCV(lithium_content[i], x_Li_max_AM, AM_material);
    } else if (type[i] == CB_type) {
      equilibrium_potential[i] = compute_OCV(lithium_content[i], x_Li_max_CB, CB_material);
    } else if (type[i] == LM_type) {
      equilibrium_potential[i] = U_OCV_LM;  // Constant for lithium metal
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixEquilibriumPotential::compute_scalar()
{
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double sum = 0.0;
  int count = 0;
  
  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit && (type[i] == AM_type || type[i] == CB_type)) {
      sum += equilibrium_potential[i];
      count++;
    }
  }
  
  double all_sum;
  int all_count;
  MPI_Allreduce(&sum,&all_sum,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&count,&all_count,1,MPI_INT,MPI_SUM,world);
  
  if (all_count > 0) return all_sum/all_count;
  return 0.0;
}