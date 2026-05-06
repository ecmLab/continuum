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
    Notes: For Symmetrical Cell Li / Li interface diffusion also considering current from SE to Li

    Input script syntax:
      fix ID group-ID lithium_diffusion keyword value ...

    Keywords:
      Li_type atype ctype       - atom types for anode and cathode Li (default 3 2)
      c_li_max value            - max Li concentration in mol/m³ (default 77101.002)
      D_Li value                - Li self-diffusion coefficient in m²/s (default 1.6e-14)
      Omega_Li value            - molar volume of Li in m³/mol (default 12.97e-6)
      s_factor value            - electrochemical time scaling factor (default 2.0e10)
      time_conv value           - time unit conversion factor to seconds (default 1.0e-6, i.e. microseconds)
      clamp_content value value - min/max bounds for lithium content (default: no clamping)
------------------------------------------------------------------------- */

#include "fix_lithium_diffusion.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "error.h"
#include "force.h"
#include "fix_property_atom.h"
#include "fix_property_atom_lithium_content.h"
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
  F(96485.0),                    // C/mol  (Faraday constant)
  c_li_max(77101.002),           // mol/m³ (Based on molar volume of Li = 0.00001297 m³/mol)
  D_Li(1.6e-14),                 // m²/s   (Li-metal self-diffusion coefficient)
  Omega_Li(12.97e-6),            // m³/mol (Molar volume of lithium)
  s_factor(2.0e10),              // Electrochemical time scaling factor (0.1s / 5e-12s)
  time_conv(1.0e-6),             // Time unit conversion to seconds (us -> s)
  clamp_content(false),          // Whether to clamp lithium content
  clamp_min(0.0),                // Min lithium content (if clamping)
  clamp_max(1.0),                // Max lithium content (if clamping)
  initial_lithium_content(0.0),
  target_lithium_content(1.0),
  max_lithium_content(1.0),
  lithium_content(NULL),
  lithium_concentration(NULL),
  current_Li_SE(NULL),
  diffusion_coefficient(NULL),
  lithium_flux(NULL),
  li_mols(NULL),
  fix_lithium_content(NULL),
  fix_lithium_concentration(NULL),
  fix_current_Li_SE(NULL),
  fix_diffusion_coefficient(NULL),
  fix_lithium_flux(NULL),
  fix_li_mols(NULL),
  fix_lithium_content_manager(NULL),
  Li_Atype(3),
  Li_Ctype(2),
  list(NULL)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix lithium_diffusion command");

  // Parse arguments
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"Li_type") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command: Li_type requires 2 values");
      Li_Atype = force->inumeric(FLERR,arg[iarg+1]);
      Li_Ctype = force->inumeric(FLERR,arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"c_li_max") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command: c_li_max requires 1 value");
      c_li_max = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"D_Li") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command: D_Li requires 1 value");
      D_Li = force->numeric(FLERR,arg[iarg+1]);
      if (D_Li <= 0.0)
        error->all(FLERR,"D_Li must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg],"Omega_Li") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command: Omega_Li requires 1 value");
      Omega_Li = force->numeric(FLERR,arg[iarg+1]);
      if (Omega_Li <= 0.0)
        error->all(FLERR,"Omega_Li must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg],"s_factor") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command: s_factor requires 1 value");
      s_factor = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"time_conv") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command: time_conv requires 1 value");
      time_conv = force->numeric(FLERR,arg[iarg+1]);
      if (time_conv <= 0.0)
        error->all(FLERR,"time_conv must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg],"clamp_content") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command: clamp_content requires 2 values (min max)");
      clamp_content = true;
      clamp_min = force->numeric(FLERR,arg[iarg+1]);
      clamp_max = force->numeric(FLERR,arg[iarg+2]);
      if (clamp_min >= clamp_max)
        error->all(FLERR,"clamp_content min must be less than max");
      iarg += 3;
    } else error->all(FLERR,"Illegal fix lithium_diffusion command: unknown keyword");
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

  // Register property/atom for lithium mols
  fix_li_mols = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumMols","property/atom","scalar",0,0,style,false));
  if(!fix_li_mols) {
    const char* fixarg[10];
    fixarg[0]="lithiumMols";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="lithiumMols";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_li_mols = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
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
    
  fix_current_Li_SE = static_cast<FixPropertyAtom*>(modify->find_fix_property("currentLiSE","property/atom","scalar",0,0,style));
  if(!fix_current_Li_SE)
    error->all(FLERR,"Fix lithium_diffusion requires currentLiSE property");
    
  fix_diffusion_coefficient = static_cast<FixPropertyAtom*>(modify->find_fix_property("diffusionCoefficient","property/atom","scalar",0,0,style));
  fix_lithium_flux = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumFlux","property/atom","scalar",0,0,style));
  fix_li_mols = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumMols","property/atom","scalar",0,0,style));
  
  if(!fix_diffusion_coefficient || !fix_lithium_flux || !fix_li_mols)
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
    if (mask[i] & groupbit && (type[i] == Li_Atype || type[i] == Li_Ctype)) {
      double radius_m = radius[i] * time_conv; // Convert radius using length scale (same as time_conv for micro units)
      li_mols[i] = ((4.0/3.0) * M_PI * radius_m * radius_m * radius_m) / Omega_Li; // mols of Li that fit in particle volume
      
      // Initialize lithium concentration
      double x_li = lithium_content[i];
      lithium_concentration[i] = x_li * c_li_max / max_lithium_content; // mol/m³
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
  current_Li_SE = fix_current_Li_SE->vector_atom;
  diffusion_coefficient = fix_diffusion_coefficient->vector_atom;
  lithium_flux = fix_lithium_flux->vector_atom;
  li_mols = fix_li_mols->vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusion::calculate_diffusion_coefficient()
{
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && (type[i] == Li_Atype || type[i] == Li_Ctype)) {
      // Use the user-specified diffusion coefficient
      // Can be extended to composition-dependent models in the future
      diffusion_coefficient[i] = D_Li; // m²/s
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
  double dt = update->dt * time_conv; // Convert timestep to seconds using user-specified conversion

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // Reset flux
  for (i = 0; i < nlocal; i++) {
    lithium_flux[i] = 0.0;
  }
  
  // Calculate diffusion flux between Li particles
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    
    if (type[i] != Li_Atype && type[i] != Li_Ctype) continue;
    
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      
      if (type[j] == Li_Atype || type[j] == Li_Ctype) {
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        r = sqrt(rsq); // simulation units (um)
        double r_SI = r * time_conv;  // Convert to m (reusing time_conv as length scale for micro units)
        
        if (r < (radius[i] + radius[j])) {
          // Calculate contact area
          double contact_area = calculate_contact_area(i, j);
          
          // Diffusion flux based on concentration gradient
          double c_diff = lithium_concentration[j] - lithium_concentration[i]; // mol/m³
          double D_eff = diffusion_coefficient[i]; // m²/s
          double flux = contact_area * D_eff * c_diff / r_SI; // mol/s
          
          lithium_flux[i] += flux; // mol/s
        }
      }
    }
    
    // Add current contribution (Equation 10)
    // Current from Li to SE contributes to lithium flux
    lithium_flux[i] -= current_Li_SE[i] / F; // mol/s (current_Li_SE is in A [C/s], F in C/mol)
  }
  
  // Update lithium content based on flux
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && (type[i] == Li_Atype || type[i] == Li_Ctype)) {
      // dn_Li/dt = flux (Equation 10)
      // s_factor scales the electrochemical time relative to the DEM timestep
      lithium_content[i] += (lithium_flux[i] * dt * s_factor) / li_mols[i]; // Li Molar Ratio
      li_mols[i] += (lithium_flux[i] * dt * s_factor); // Li Mols in Li Metal
 
      // Clamp lithium content if user requested
      if (clamp_content) {
        if (lithium_content[i] < clamp_min) lithium_content[i] = clamp_min;
        if (lithium_content[i] > clamp_max) lithium_content[i] = clamp_max;
      }
      
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
  double length_conversion = time_conv;  // Reuse time_conv as length scale for micro units
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
    if (type[i] == Li_Atype || type[i] == Li_Ctype) {
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