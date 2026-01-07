/* ----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine
    http://www.dcs-computing.com, office@dcs-computing.com
    LIGGGHTS® is part of CFDEM®project: http://www.liggghts.com | http://www.cfdem.com

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    Created: Joseph Vazquez Mercado, RIT 2025
    Copyright 2024-     DCS Computing GmbH, Linz
    Notes: Radius expansion for metal-carbon lithium alloying system
           V = Omega_metal * n_metal + Omega_Li_eff * n_Li
           where:
           - Omega_metal = 10.96e-6 m³/mol (molar volume of metal, e.g., Ag)
           - Omega_Li_eff = 9.0e-6 m³/mol (effective molar volume of Li in alloy)
------------------------------------------------------------------------- */

#include "fix_metal_radius_expansion.h"
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

FixMetalRadiusExpansion::FixMetalRadiusExpansion(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  Omega_metal(10.96e-6),     // m³/mol molar volume of metal (e.g., Ag)
  Omega_Li_eff(9.0e-6),      // m³/mol effective molar volume of Li in metal alloy
  update_frequency(10),
  update_count(0),
  next_update_step(0),
  manual_trigger(false),
  lithium_content(NULL),
  li_mols(NULL),
  metal_mols(NULL),
  initial_radius(NULL),
  particle_volume(NULL),
  fix_lithium_content(NULL),
  fix_li_mols(NULL),
  fix_metal_mols(NULL),
  fix_initial_radius(NULL),
  fix_particle_volume(NULL),
  metal_type(1),
  carbon_type(2)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix metal_radius_expansion command");

  // Parse arguments
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"metal_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix metal_radius_expansion command");
      metal_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"carbon_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix metal_radius_expansion command");
      carbon_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"update_every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix metal_radius_expansion command");
      update_frequency = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"Omega_metal") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix metal_radius_expansion command");
      Omega_metal = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"Omega_Li_eff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix metal_radius_expansion command");
      Omega_Li_eff = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix metal_radius_expansion command");
  }

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;

  vector_flag = 1;
  size_vector = 3;
  extvector = 0;
}

/* ---------------------------------------------------------------------- */

void FixMetalRadiusExpansion::post_create()
{
  // Register property/atom for initial radius
  fix_initial_radius = static_cast<FixPropertyAtom*>(modify->find_fix_property("initialRadius","property/atom","scalar",0,0,style,false));
  if(!fix_initial_radius) {
    const char* fixarg[9];
    fixarg[0]="initialRadius";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="initialRadius";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_initial_radius = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  // Register property/atom for particle volume
  fix_particle_volume = static_cast<FixPropertyAtom*>(modify->find_fix_property("particleVolume","property/atom","scalar",0,0,style,false));
  if(!fix_particle_volume) {
    const char* fixarg[9];
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

FixMetalRadiusExpansion::~FixMetalRadiusExpansion()
{
}

/* ---------------------------------------------------------------------- */

int FixMetalRadiusExpansion::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMetalRadiusExpansion::init()
{
  // Find required fixes
  fix_lithium_content = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumContent","property/atom","scalar",0,0,style));
  if(!fix_lithium_content)
    error->all(FLERR,"Fix metal_radius_expansion requires lithiumContent property");
  
  fix_li_mols = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumMols","property/atom","scalar",0,0,style));
  if(!fix_li_mols)
    error->all(FLERR,"Fix metal_radius_expansion requires lithiumMols property");

  fix_metal_mols = static_cast<FixPropertyAtom*>(modify->find_fix_property("metalMols","property/atom","scalar",0,0,style));
  if(!fix_metal_mols)
    error->all(FLERR,"Fix metal_radius_expansion requires metalMols property");
    
  fix_initial_radius = static_cast<FixPropertyAtom*>(modify->find_fix_property("initialRadius","property/atom","scalar",0,0,style));
  fix_particle_volume = static_cast<FixPropertyAtom*>(modify->find_fix_property("particleVolume","property/atom","scalar",0,0,style));
  
  if(!fix_initial_radius || !fix_particle_volume)
    error->all(FLERR,"Could not find required property/atom fixes");

  updatePtrs();
  
  // Set next update step
  next_update_step = update->ntimestep + update_frequency;
}

/* ---------------------------------------------------------------------- */

void FixMetalRadiusExpansion::setup(int vflag)
{
  updatePtrs();
  
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      // Store initial radius (units: μm)
      initial_radius[i] = radius[i];
      
      // Calculate initial particle volume (m³)
      double radius_m = radius[i] * 1.0e-6;  // Convert μm to m
      double V_initial = (4.0/3.0) * M_PI * radius_m * radius_m * radius_m;
      particle_volume[i] = V_initial;
      
      if (type[i] == metal_type) {
        // For metal particles: initially n_Li = 0
        // V_initial = Omega_metal * n_metal
        // Therefore: n_metal = V_initial / Omega_metal
        metal_mols[i] = V_initial / Omega_metal;
        
        // Verify li_mols is initialized to 0 (should be done by diffusion fix)
        if (li_mols[i] < 0.0) li_mols[i] = 0.0;
      }
      // Carbon particles don't expand - they maintain constant size
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMetalRadiusExpansion::end_of_step()
{
  // No max_updates limit - always update when frequency is reached
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

void FixMetalRadiusExpansion::updatePtrs()
{
  lithium_content = fix_lithium_content->vector_atom;
  li_mols = fix_li_mols->vector_atom;
  metal_mols = fix_metal_mols->vector_atom;
  initial_radius = fix_initial_radius->vector_atom;
  particle_volume = fix_particle_volume->vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixMetalRadiusExpansion::update_particle_radius()
{
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && type[i] == metal_type) {
      double n_metal = metal_mols[i];  // Moles of metal (constant)
      double n_Li = li_mols[i];        // Moles of lithium (increases with charging)
      
      // Volume equation: V = Omega_metal * n_metal + Omega_Li_eff * n_Li
      particle_volume[i] = Omega_metal * n_metal + Omega_Li_eff * n_Li;  // m³
      
      // Calculate new radius from volume
      // V = (4/3) * pi * r³  =>  r = cbrt(3V / (4*pi))
      double new_radius_m = cbrt((3.0 * particle_volume[i]) / (4.0 * M_PI));  // m
      double new_radius_um = new_radius_m * 1.0e6;  // Convert to μm
      
      // No radius cap - allow unlimited expansion
      // Only ensure radius is positive and physically reasonable
      if (new_radius_um > 0.0) {
        radius[i] = new_radius_um;
      }
    }
    // Carbon particles don't change radius - they are just pathways
  }
  
  // Force complete neighbor list rebuild
  neighbor->delay = 0;
  neighbor->every = 1;
  neighbor->ago = neighbor->every;

  // Forward communication to update ghost particles
  comm->forward_comm();
}

/* ---------------------------------------------------------------------- */

void FixMetalRadiusExpansion::trigger_radius_update()
{
  manual_trigger = true;
}

/* ---------------------------------------------------------------------- */

double FixMetalRadiusExpansion::compute_scalar()
{
  // Return current update count
  return static_cast<double>(update_count);
}

/* ---------------------------------------------------------------------- */

double FixMetalRadiusExpansion::compute_vector(int n)
{
  if (n == 0) {
    // Return update count
    return static_cast<double>(update_count);
  } else if (n == 1) {
    // Return average expansion ratio for metal particles
    int nlocal = atom->nlocal;
    int *type = atom->type;
    int *mask = atom->mask;
    double *radius = atom->radius;
    double sum_ratio = 0.0;
    int count = 0;
    
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit && type[i] == metal_type && initial_radius[i] > 0.0) {
        sum_ratio += radius[i] / initial_radius[i];
        count++;
      }
    }
    
    double all_sum, all_count_d;
    int all_count;
    MPI_Allreduce(&sum_ratio, &all_sum, 1, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(&count, &all_count, 1, MPI_INT, MPI_SUM, world);
    
    if (all_count > 0) return all_sum / all_count;
    return 1.0;
  } else if (n == 2) {
    // Return average Li/Metal ratio
    int nlocal = atom->nlocal;
    int *type = atom->type;
    int *mask = atom->mask;
    double sum_ratio = 0.0;
    int count = 0;
    
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit && type[i] == metal_type) {
        sum_ratio += lithium_content[i];
        count++;
      }
    }
    
    double all_sum;
    int all_count;
    MPI_Allreduce(&sum_ratio, &all_sum, 1, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(&count, &all_count, 1, MPI_INT, MPI_SUM, world);
    
    if (all_count > 0) return all_sum / all_count;
    return 0.0;
  }
  return 0.0;
}