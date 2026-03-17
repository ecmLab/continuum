/* ----------------------------------------------------------------------
    LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
    Transfer Simulations

    Fix for managing lithium content in Si particles for battery simulation

    Created: Joseph Vazquez Mercado, RIT 2025
------------------------------------------------------------------------- */

#include "fix_property_atom_lithium_content.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "error.h"
#include "force.h"
#include "fix_property_atom.h"
#include <cstring>
#include "comm.h"
#include "memory.h"
#include <cstdio>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

void FixPropertyAtomLithiumContent::post_create()
{
  // Register property/atom for lithium content
  fix_lithium_content = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumContent","property/atom","scalar",0,0,style,false));
  if(!fix_lithium_content) {
    const char* fixarg[10];
    fixarg[0]="lithiumContent";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="lithiumContent";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    char arg8[30];
    sprintf(arg8,"%f",initial_lithium_content);
    fixarg[8]=arg8;
    fix_lithium_content = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }
}

/* ---------------------------------------------------------------------- */

FixPropertyAtomLithiumContent::FixPropertyAtomLithiumContent(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  initial_lithium_content(0.0),
  target_lithium_content(0.9),
  max_lithium_content(1.0),
  lithium_content(NULL),
  fix_lithium_content(NULL)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix property/atom/lithium_content command");

  // Parse optional arguments
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"initial") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix property/atom/lithium_content command");
      initial_lithium_content = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"target") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix property/atom/lithium_content command");
      target_lithium_content = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"max") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix property/atom/lithium_content command");
      max_lithium_content = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix property/atom/lithium_content command");
  }

  // Validate input
  if (initial_lithium_content < 0.0 || initial_lithium_content > max_lithium_content)
    error->all(FLERR,"Invalid initial lithium content");
  if (target_lithium_content < 0.0 || target_lithium_content > max_lithium_content)
    error->all(FLERR,"Invalid target lithium content");
  if (max_lithium_content <= 0.0)
    error->all(FLERR,"Invalid maximum lithium content");

  // Don't create property/atom in constructor - wait for post_create()
  // This avoids issues with modify->add_fix() during construction

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;

  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 0;
}

/* ---------------------------------------------------------------------- */

FixPropertyAtomLithiumContent::~FixPropertyAtomLithiumContent()
{
}

/* ---------------------------------------------------------------------- */

int FixPropertyAtomLithiumContent::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtomLithiumContent::init()
{
  // Find the lithium content property
  fix_lithium_content = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumContent","property/atom","scalar",0,0,style));
  if(!fix_lithium_content)
    error->one(FLERR,"Could not find fix property/atom lithiumContent");
    
  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtomLithiumContent::setup(int vflag)
{
  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtomLithiumContent::min_setup(int vflag)
{
  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtomLithiumContent::updatePtrs()
{
  lithium_content = fix_lithium_content->vector_atom;
}

/* ---------------------------------------------------------------------- */

double FixPropertyAtomLithiumContent::compute_scalar()
{
  // Return average lithium content
  int nlocal = atom->nlocal;
  double sum = 0.0;
  int count = 0;
  
  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      sum += lithium_content[i];
      count++;
    }
  }
  
  double all_sum, all_count;
  MPI_Allreduce(&sum,&all_sum,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&count,&all_count,1,MPI_INT,MPI_SUM,world);
  
  if (all_count > 0) return all_sum/all_count;
  return 0.0;
}

/* ---------------------------------------------------------------------- */

double FixPropertyAtomLithiumContent::compute_vector(int n)
{
  if (n == 0) return initial_lithium_content;
  else if (n == 1) return target_lithium_content;
  else if (n == 2) return max_lithium_content;
  return 0.0;
}
