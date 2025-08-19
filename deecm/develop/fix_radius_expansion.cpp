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
  Omega_Si(10.96e-6),        // m³/mol
  Omega_Li_eff(9e-6),        // m³/mol
  update_frequency(10),
  max_updates(10),
  update_count(0),
  next_update_step(0),
  manual_trigger(false),
  lithium_content(NULL),
  silicon_content(NULL),
  initial_radius(NULL),
  particle_volume(NULL),
  fix_lithium_content(NULL),
  fix_silicon_content(NULL),
  fix_initial_radius(NULL),
  fix_particle_volume(NULL),
  AM_type(1)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix radius_expansion command");

  // Parse arguments
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"AM_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix radius_expansion command");
      AM_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"update_every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix radius_expansion command");
      update_frequency = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"max_updates") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix radius_expansion command");
      max_updates = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix radius_expansion command");
  }

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;

  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 0;
}

/* ---------------------------------------------------------------------- */

void FixRadiusExpansion::post_create()
{
  // Register property/atom for initial radius (to be compatible with lithium_diffusion)
  fix_initial_radius = static_cast<FixPropertyAtom*>(modify->find_fix_property("initialRadius","property/atom","scalar",0,0,style,false));
  if(!fix_initial_radius) {
    const char* fixarg[10];
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
  // Find required fixes
  fix_lithium_content = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumContent","property/atom","scalar",0,0,style));
  if(!fix_lithium_content)
    error->all(FLERR,"Fix radius_expansion requires lithiumContent property");
  
  // Find silicon content property (created by lithium_diffusion)
  fix_silicon_content = static_cast<FixPropertyAtom*>(modify->find_fix_property("siliconContent","property/atom","scalar",0,0,style));
  if(!fix_silicon_content)
    error->all(FLERR,"Fix radius_expansion requires siliconContent property (should be created by fix lithium_diffusion)");
    
  fix_initial_radius = static_cast<FixPropertyAtom*>(modify->find_fix_property("initialRadius","property/atom","scalar",0,0,style));
  fix_particle_volume = static_cast<FixPropertyAtom*>(modify->find_fix_property("particleVolume","property/atom","scalar",0,0,style));
  
  if(!fix_initial_radius || !fix_particle_volume)
    error->all(FLERR,"Could not find required property/atom fixes");

  updatePtrs();
  
  // Set next update step
  next_update_step = update->ntimestep + update_frequency;
}

/* ---------------------------------------------------------------------- */

void FixRadiusExpansion::setup(int vflag)
{
  updatePtrs();
  
  // Store initial radius values that were set by lithium_diffusion
  // lithium_diffusion should have already run its setup() and populated silicon_content
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  // Store the initial radius (this is the radius at simulation start)
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && type[i] == AM_type) {
      // Only set initial_radius if it hasn't been set yet (value is 0.0)
      if (initial_radius[i] == 0.0) {
        initial_radius[i] = radius[i];
      }
      
      // Calculate initial particle volume
      double V_initial = (4.0/3.0) * M_PI * initial_radius[i] * initial_radius[i] * initial_radius[i];
      particle_volume[i] = V_initial;
    }
  }
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
  silicon_content = fix_silicon_content->vector_atom;
  initial_radius = fix_initial_radius->vector_atom;
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
    if (mask[i] & groupbit && type[i] == AM_type) {
      double x_Li = lithium_content[i];  // Li/Si molar ratio
      double n_Si = silicon_content[i];  // Moles of silicon
      double n_Li = x_Li * n_Si;          // Moles of lithium

      // Final volume of the particle: V = Omega_Si * n_Si + Omega_Li * n_Li
      particle_volume[i] = Omega_Si * n_Si + Omega_Li_eff * n_Li;

      // The new radius is derived from the new volume
      // V = (4/3) * pi * r^3, so r = cbrt(3V/(4*pi))
      radius[i] = cbrt((3.0 * particle_volume[i]) / (4.0 * M_PI));

      // Convert from m to μm (LIGGGHTS units)
      radius[i] = radius[i] * 1.0e6;

      // Safety check - don't let radius shrink below initial
      if (radius[i] < initial_radius[i]) {
        radius[i] = initial_radius[i];
      } else if (radius[i] > 2.0 * initial_radius[i]) {
        radius[i] = 2.0 * initial_radius[i];
        if (comm->me == 0 && screen) {
          fprintf(screen,"Warning: Particle %d reached maximum expansion limit\n", atom->tag[i]);
        }
      }
    }  
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
  // Return current update count
  return (double)update_count;
}

/* ---------------------------------------------------------------------- */

double FixRadiusExpansion::compute_vector(int n)
{
  if (n == 0) return (double)update_count;
  else if (n == 1) return (double)(max_updates - update_count);
  return 0.0;
}
