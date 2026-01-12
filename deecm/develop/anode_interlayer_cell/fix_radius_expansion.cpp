/* ----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine
    Contributing author: Joseph Vazquez Mercado, RIT 2025
    Copyright 2024-     DCS Computing GmbH, Linz
    Notes: Radius expansion for AM and CB particles based on Li content
------------------------------------------------------------------------- */

#include "fix_radius_expansion.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "error.h"
#include "force.h"
#include "fix_property_atom.h"
#include "neighbor.h"
#include "comm.h"
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixRadiusExpansion::FixRadiusExpansion(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  Omega_Li_AM(6.0593e-6),         // m³/mol (effective molar volume of AgLi4 - Based on max concentration)
  Omega_Li_CB(4.9307e-6),         // m³/mol (effective molar volume in CB - Based on max concentration)
  update_frequency(10),
  max_updates(10),
  update_count(0),
  next_update_step(0),
  manual_trigger(false),
  lithium_content(NULL),
  li_mols(NULL),
  initial_volume(NULL),         // Changed from initial_radius
  particle_volume(NULL),
  fix_lithium_content(NULL),
  fix_li_mols(NULL),
  fix_initial_volume(NULL),     // Changed from fix_initial_radius
  fix_particle_volume(NULL),
  AM_type(1),
  CB_type(2)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix radius_expansion command");

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"AM_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix radius_expansion command");
      AM_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"CB_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix radius_expansion command");
      CB_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"update_every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix radius_expansion command");
      update_frequency = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"max_updates") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix radius_expansion command");
      max_updates = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"Omega_Li_AM") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix radius_expansion command");
      Omega_Li_AM = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"Omega_Li_CB") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix radius_expansion command");
      Omega_Li_CB = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix radius_expansion command");
  }

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;

  vector_flag = 1;
  size_vector = 2;
  extvector = 0;
}

/* ---------------------------------------------------------------------- */

void FixRadiusExpansion::post_create()
{
  // Register property/atom for particle volume
  fix_particle_volume = static_cast<FixPropertyAtom*>(modify->find_fix_property("particleVolume","property/atom","scalar",0,0,style,false));
  if(!fix_particle_volume) {
    const char* fixarg[10];
    fixarg[0]="particleVolume";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="particleVolume";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_particle_volume = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }
}

/* ---------------------------------------------------------------------- */

FixRadiusExpansion::~FixRadiusExpansion()
{
}

/* ---------------------------------------------------------------------- */

int FixRadiusExpansion::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRadiusExpansion::init()
{
  fix_lithium_content = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumContent","property/atom","scalar",0,0,style));
  if(!fix_lithium_content)
    error->all(FLERR,"Fix radius_expansion requires lithiumContent property");
  
  fix_li_mols = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumMols","property/atom","scalar",0,0,style));
  if(!fix_li_mols)
    error->all(FLERR,"Fix radius_expansion requires lithiumMols property (created by fix lithium_diffusion)");
  
  // Find initialVolume (managed by fix lithium_diffusion)
  fix_initial_volume = static_cast<FixPropertyAtom*>(modify->find_fix_property("initialVolume","property/atom","scalar",0,0,style));
  if(!fix_initial_volume)
    error->all(FLERR,"Fix radius_expansion requires initialVolume property (created by fix lithium_diffusion)"); 
    
  fix_particle_volume = static_cast<FixPropertyAtom*>(modify->find_fix_property("particleVolume","property/atom","scalar",0,0,style));
  if(!fix_particle_volume)
    error->all(FLERR,"Could not find particleVolume property/atom fix");

  updatePtrs();
  
  next_update_step = update->ntimestep + update_frequency;
}

/* ---------------------------------------------------------------------- */

void FixRadiusExpansion::setup(int vflag)
{
  updatePtrs();
  
  // Calculate the initial state based on initialVolume and whatever 
  // lithium mols exist at step 0 (from fix lithium_diffusion)
  update_particle_radius();
}

/* ---------------------------------------------------------------------- */

void FixRadiusExpansion::end_of_step()
{
  if (update_count >= max_updates) return;
  
  if (manual_trigger || update->ntimestep >= next_update_step) {
    update_particle_radius();
    update_count++;
    next_update_step = update->ntimestep + update_frequency;
    manual_trigger = false;
    
    // Trigger neighbor list rebuild
    neighbor->delay = 0;
    neighbor->ago = neighbor->every;
  }
}

/* ---------------------------------------------------------------------- */

void FixRadiusExpansion::updatePtrs()
{
  lithium_content = fix_lithium_content->vector_atom;
  li_mols = fix_li_mols->vector_atom;
  initial_volume = fix_initial_volume->vector_atom; // Pointer to initialVolume
  particle_volume = fix_particle_volume->vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixRadiusExpansion::update_particle_radius()
{
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    
    int type_i = type[i];
    
    // Only update AM and CB particles
    if (type_i != AM_type && type_i != CB_type) continue;
    
    double n_Li = li_mols[i];  // Moles of lithium from diffusion fix
    double Omega_eff;
    
    // Select appropriate molar volume based on particle type
    if (type_i == AM_type) {
      Omega_eff = Omega_Li_AM;
    } else {
      Omega_eff = Omega_Li_CB;
    }
    
    // Calculate new volume based on lithium content
    // V = initial_volume[i] + Omega_Li_eff * n_Li (in m³)
    particle_volume[i] = initial_volume[i] + Omega_eff * n_Li;
    
    // Calculate new radius from volume
    // V = (4/3) * pi * r³  =>  r = cbrt(3V / (4*pi))
    double r_m = cbrt((3.0 * particle_volume[i]) / (4.0 * M_PI)); // in meters
    
    // Convert from m to μm (LIGGGHTS units)
    radius[i] = r_m * 1.0e6;
  }
  
  // Force complete neighbor list rebuild
  neighbor->delay = 0;
  neighbor->every = 1;
  neighbor->ago = neighbor->every;

  // Forward communication to update ghost particles
  comm->forward_comm();
}

/* ---------------------------------------------------------------------- */

void FixRadiusExpansion::trigger_radius_update()
{
  manual_trigger = true;
}

/* ---------------------------------------------------------------------- */

double FixRadiusExpansion::compute_scalar()
{
  return (double)update_count;
}

/* ---------------------------------------------------------------------- */

double FixRadiusExpansion::compute_vector(int n)
{
  if (n == 0) return (double)update_count;
  else if (n == 1) return (double)(max_updates - update_count);
  return 0.0;
}