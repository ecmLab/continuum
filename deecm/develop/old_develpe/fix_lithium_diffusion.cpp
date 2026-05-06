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

#include "fix_lithium_diffusion.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "error.h"
#include "force.h"
#include "fix_property_atom.h"
#include "fix_property_atom_lithium_content.h"
#include "fix_exchange_current_density.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixLithiumDiffusion::FixLithiumDiffusion(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  D_min(1e-15),              // m²/s
  D_poor(1e-14),             // m²/s  
  D_rich(1e-13),             // m²/s
  kappa(20),
  F(96485.0),                // C/mol
  c_li_max(83874.0),         // mol/m³
  initial_lithium_content(0.0),  // Default value
  target_lithium_content(0.9),    // Default value
  max_lithium_content(1.0),       // Default value
  lithium_content(NULL),
  lithium_concentration(NULL),
  current_AM_SE(NULL),
  diffusion_coefficient(NULL),
  lithium_flux(NULL),
  silicon_content(NULL),
  fix_lithium_content(NULL),
  fix_lithium_concentration(NULL),
  fix_current_AM_SE(NULL),
  fix_diffusion_coefficient(NULL),
  fix_lithium_flux(NULL),
  fix_silicon_content(NULL),
  fix_lithium_content_manager(NULL),
  AM_type(1),
  list(NULL)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix lithium_diffusion command");

  // Parse arguments
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"AM_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      AM_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"D_min") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      D_min = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"D_poor") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      D_poor = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"D_rich") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      D_rich = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"kappa") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      kappa = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"c_li_max") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      c_li_max = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix lithium_diffusion command");
  }

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusion::post_create()
{
  // Register property/atom for diffusion coefficient
  fix_diffusion_coefficient = static_cast<FixPropertyAtom*>(modify->find_fix_property("diffusionCoefficient","property/atom","scalar",0,0,style,false));
  if(!fix_diffusion_coefficient) {
    const char* fixarg[10];
    fixarg[0]="diffusionCoefficient";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="diffusionCoefficient";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_diffusion_coefficient = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  // Register property/atom for lithium flux
  fix_lithium_flux = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumFlux","property/atom","scalar",0,0,style,false));
  if(!fix_lithium_flux) {
    const char* fixarg[10];
    fixarg[0]="lithiumFlux";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="lithiumFlux";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_lithium_flux = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  // Register property/atom for silicon content
  fix_silicon_content = static_cast<FixPropertyAtom*>(modify->find_fix_property("siliconContent","property/atom","scalar",0,0,style,false));
  if(!fix_silicon_content) {
    const char* fixarg[10];
    fixarg[0]="siliconContent";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="siliconContent";
    fixarg[4]="scalar";
    fixarg[5]="no";    // no restart
    fixarg[6]="yes";   // communicate
    fixarg[7]="no";    // no borders
    fixarg[8]="0.0";
    fix_silicon_content = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }
}

/* ---------------------------------------------------------------------- */

FixLithiumDiffusion::~FixLithiumDiffusion()
{
}

/* ---------------------------------------------------------------------- */

int FixLithiumDiffusion::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusion::init()
{
  // Find required fixes
  fix_lithium_content = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumContent","property/atom","scalar",0,0,style));
  if(!fix_lithium_content)
    error->all(FLERR,"Fix lithium_diffusion requires lithiumContent property");
    
  fix_lithium_concentration = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumConcentration","property/atom","scalar",0,0,style));
  if(!fix_lithium_concentration)
    error->all(FLERR,"Fix lithium_diffusion requires lithiumConcentration property");
    
  fix_current_AM_SE = static_cast<FixPropertyAtom*>(modify->find_fix_property("currentAMSE","property/atom","scalar",0,0,style));
  if(!fix_current_AM_SE)
    error->all(FLERR,"Fix lithium_diffusion requires currentAMSE property");
    
  fix_diffusion_coefficient = static_cast<FixPropertyAtom*>(modify->find_fix_property("diffusionCoefficient","property/atom","scalar",0,0,style));
  fix_lithium_flux = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumFlux","property/atom","scalar",0,0,style));
  fix_silicon_content = static_cast<FixPropertyAtom*>(modify->find_fix_property("siliconContent","property/atom","scalar",0,0,style));
  
  if(!fix_diffusion_coefficient || !fix_lithium_flux || !fix_silicon_content)
    error->all(FLERR,"Could not find required property/atom fixes");

  // Find the lithium content manager fix to get initial/target/max values
  fix_lithium_content_manager = NULL;
  for (int i = 0; i < modify->nfix; i++) {
    if (strstr(modify->fix[i]->style,"property/atom/lithium_content")) {
      fix_lithium_content_manager = static_cast<FixPropertyAtomLithiumContent*>(modify->fix[i]);
      break;
    }
  }
  
  if(!fix_lithium_content_manager)
    error->all(FLERR,"Fix lithium_diffusion requires fix property/atom/lithium_content to be defined");
  
  // Get lithium content parameters from the manager fix
  initial_lithium_content = fix_lithium_content_manager->compute_vector(0);
  target_lithium_content = fix_lithium_content_manager->compute_vector(1);
  max_lithium_content = fix_lithium_content_manager->compute_vector(2);

  // Request neighbor list
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusion::setup(int vflag)
{
  updatePtrs();
  // Calculate initial diffusion coefficient
  calculate_diffusion_coefficient();
  
  // Store initial radius values (only done once at setup)
  double *radius = atom->radius;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  int *mask = atom->mask;
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && type[i] == AM_type) {
      double Omega_NMC = 20.44e-6; // m³/mol (Molar Volume of Silicon 10.96e-6)
      double radius_m = radius[i] * 1.0e-6;
      silicon_content[i] = ((4.0/3.0) * M_PI * radius_m * radius_m * radius_m) / ( Omega_NMC);
      
      // Initialize lithium concentration
      double x_li = lithium_content[i];
      lithium_concentration[i] = x_li * c_li_max / max_lithium_content;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusion::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusion::pre_force(int vflag)
{
  updatePtrs();
  
  // Calculate diffusion coefficient based on lithium content
  calculate_diffusion_coefficient();
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusion::post_force(int vflag)
{
  // Update lithium content based on diffusion and current
  update_lithium_content();
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusion::updatePtrs()
{
  lithium_content = fix_lithium_content->vector_atom;
  lithium_concentration = fix_lithium_concentration->vector_atom;
  current_AM_SE = fix_current_AM_SE->vector_atom;
  diffusion_coefficient = fix_diffusion_coefficient->vector_atom;
  lithium_flux = fix_lithium_flux->vector_atom;
  silicon_content = fix_silicon_content->vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusion::calculate_diffusion_coefficient()
{
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && type[i] == AM_type) {
      double x_li = lithium_content[i];
      // double x_ratio = x_li / max_lithium_content; // Chi ratio using max_lithium_content
      
      // // Equation 11: D_Li-Si = D_min + D_poor * exp(-kappa*x_Li/x_max) + D_rich * exp(kappa*(x_Li/x_max - 1))
      // double term1 = D_poor * exp(-kappa * x_ratio);
      // double term2 = D_rich * exp(kappa * (x_ratio - 1.0));
      
      // diffusion_coefficient[i] = D_min + term1 + term2;
      diffusion_coefficient[i] = (2.242731e-14 * x_li) + 2.531562e-15; // Linear increase Diffusion Coefficient m²/s based on x_li from 0 to 1 (based on DOI: 10.1021/acs.chemmater.2c03834)
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusion::update_lithium_content()
{
  int i,j,ii,jj,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,r;
  
  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt = update->dt * 1.0e-6; // us to s  

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // Reset flux
  for (i = 0; i < nlocal; i++) {
    lithium_flux[i] = 0.0;
  }
  
  // Calculate diffusion flux between AM particles
  for (ii = 0; ii < inum; ii++) { // 1
    i = ilist[ii];
    
    if (type[i] != AM_type) continue;
    
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    for (jj = 0; jj < jnum; jj++) { // 2
      j = jlist[jj];
      j &= NEIGHMASK;
      
      if (type[j] == AM_type) {  // 3
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        r = sqrt(rsq); // um
        double r_SI = r * 1.0e-6;  // Convert to m (Distance between AM_q-AM_p)
        
        if (r < (radius[i] + radius[j])) { // 4
          // Calculate contact area
          double contact_area = calculate_contact_area(i, j);
          
          // Diffusion flux based on concentration gradient
          double c_diff = lithium_concentration[j] - lithium_concentration[i];
          double D_Li_AM = diffusion_coefficient[i];
          double flux = contact_area * D_Li_AM * c_diff / r_SI;
          
          lithium_flux[i] += flux;

        } // 4
      } // 3
    } // 2 end of j
    
    // Add current contribution (Equation 10)
    // Current from AM to SE contributes to lithium flux
    lithium_flux[i] -= current_AM_SE[i] / F;
  } // 1
  
  // Update lithium content based on flux
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && type[i] == AM_type) {
      // dn_Li/dt = flux (Equation 10)
      
      // Update Notes: The scalling factor needs to be toned down initially to avoid large jumps in lithium content
      // When Lithium content is at least 0.7 then you can use the max scaling factor of 2e16
      // The scaling factor is based on a timestep of 5e-6 us in main, so it simulates 0.1s of EC in 1 run
      double s_factor = 2.0e11; // 1s / 5e-12s 
      lithium_content[i] += (lithium_flux[i] * dt * s_factor) / silicon_content[i]; // Li Molar Ratio (Silicon Content is actually NMC811)
 
      // Ensure lithium content stays within bounds using values from lithium_content manager
      if (lithium_content[i] < initial_lithium_content) lithium_content[i] = initial_lithium_content;
      if (lithium_content[i] > target_lithium_content) lithium_content[i] = target_lithium_content;
      
      // Update lithium concentration based on new lithium content
      lithium_concentration[i] = lithium_content[i] * c_li_max / max_lithium_content;
    }
  }
  
  // Forward communication
  fix_lithium_content->do_forward_comm();
  fix_lithium_concentration->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

double FixLithiumDiffusion::calculate_contact_area(int i, int j)
{
  double *radius = atom->radius;
  double **x = atom->x;
  double length_conversion = 1.0e-6;  // LAMMPS units to m
  double delx = (x[i][0] - x[j][0]) * length_conversion;
  double dely = (x[i][1] - x[j][1]) * length_conversion;
  double delz = (x[i][2] - x[j][2]) * length_conversion;
  double r = sqrt(delx*delx + dely*dely + delz*delz); // m
  double radi = radius[i] * length_conversion; // m
  double radj = radius[j] * length_conversion; // m

  double radsum = radi + radj; // m
  
  if (r >= radsum) return 0.0;
  
  // Simple overlap-based contact area
  double delta_n = radsum - r;
  double reff = radi * radj / radsum;
  return M_PI * delta_n * reff; // returns area in m2
}

/* ---------------------------------------------------------------------- */

double FixLithiumDiffusion::compute_scalar()
{
  // Return average lithium flux
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double sum = 0.0;
  int count = 0;
  
  for (int i = 0; i < nlocal; i++) {
    if (type[i] == AM_type) {
      sum += fabs(lithium_flux[i]);
      count++;
    }
  }
  
  double all_sum, all_count;
  MPI_Allreduce(&sum,&all_sum,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&count,&all_count,1,MPI_INT,MPI_SUM,world);
  
  if (all_count > 0) return all_sum/all_count;
  return 0.0;
}
