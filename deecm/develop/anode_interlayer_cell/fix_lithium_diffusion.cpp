/* ----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine
    Contributing author: Joseph Vazquez Mercado, RIT 2025
    Copyright 2024-     DCS Computing GmbH, Linz
    Notes: Multi-material Li diffusion: AM, CB, and LM (constant source)
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
  F(96485.0),                    // C/mol
  c_li_max(77101.002),           // mol/m³ (LM constant concentration)
  Omega_Li_LM(13.02e-6),         // m³/mol (Li molar volume in LM)
  Omega_Li_AM(8.91e-6),          // m³/mol (Li molar volume in AM)
  Omega_Li_CB(8.91e-6),          // m³/mol (Li molar volume in CB)
  initial_lithium_content(0.0),
  target_lithium_content(1.0),
  max_lithium_content(1.0),
  // Diffusion coefficients (m²/s)
  D_AM_AM(1.0e-15),
  D_AM_CB(1.0e-15),
  D_AM_LM(1.0e-15),
  D_CB_AM(1.0e-15),
  D_CB_CB(1.0e-15),
  D_CB_LM(1.0e-15),
  lithium_content(NULL),
  lithium_concentration(NULL),
  current_SE_Li(NULL),
  diffusion_coefficient(NULL),
  lithium_flux(NULL),
  li_mols(NULL),
  fix_lithium_content(NULL),
  fix_lithium_concentration(NULL),
  fix_current_SE_Li(NULL),
  fix_diffusion_coefficient(NULL),
  fix_lithium_flux(NULL),
  fix_li_mols(NULL),
  fix_lithium_content_manager(NULL),
  AM_type(1),
  CB_type(2),
  LM_type(4),
  list(NULL)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix lithium_diffusion command");

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"AM_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      AM_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"CB_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      CB_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"LM_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      LM_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"c_li_max") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      c_li_max = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"D_AM_AM") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      D_AM_AM = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"D_AM_CB") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      D_AM_CB = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"D_AM_LM") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      D_AM_LM = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"D_CB_AM") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      D_CB_AM = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"D_CB_CB") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      D_CB_CB = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"D_CB_LM") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      D_CB_LM = force->numeric(FLERR,arg[iarg+1]);
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
  fix_lithium_content = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumContent","property/atom","scalar",0,0,style));
  if(!fix_lithium_content)
    error->all(FLERR,"Fix lithium_diffusion requires lithiumContent property");
    
  fix_lithium_concentration = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumConcentration","property/atom","scalar",0,0,style));
  if(!fix_lithium_concentration)
    error->all(FLERR,"Fix lithium_diffusion requires lithiumConcentration property");
    
  fix_current_SE_Li = static_cast<FixPropertyAtom*>(modify->find_fix_property("currentSELi","property/atom","scalar",0,0,style,false));
  // current_SE_Li is optional - may not exist for pure diffusion simulations
    
  fix_diffusion_coefficient = static_cast<FixPropertyAtom*>(modify->find_fix_property("diffusionCoefficient","property/atom","scalar",0,0,style));
  fix_lithium_flux = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumFlux","property/atom","scalar",0,0,style));
  fix_li_mols = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumMols","property/atom","scalar",0,0,style));
  
  if(!fix_diffusion_coefficient || !fix_lithium_flux || !fix_li_mols)
    error->all(FLERR,"Could not find required property/atom fixes");

  fix_lithium_content_manager = NULL;
  for (int i = 0; i < modify->nfix; i++) {
    if (strstr(modify->fix[i]->style,"property/atom/lithium_content")) {
      fix_lithium_content_manager = static_cast<FixPropertyAtomLithiumContent*>(modify->fix[i]);
      break;
    }
  }
  
  if(fix_lithium_content_manager) {
    initial_lithium_content = fix_lithium_content_manager->compute_vector(0);
    target_lithium_content = fix_lithium_content_manager->compute_vector(1);
    max_lithium_content = fix_lithium_content_manager->compute_vector(2);
  }

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
  
  double *radius = atom->radius;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  int *mask = atom->mask;
  
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    
    double radius_m = radius[i] * 1.0e-6;
    double volume = (4.0/3.0) * M_PI * radius_m * radius_m * radius_m;
    
    if (type[i] == AM_type) {
      li_mols[i] = volume / Omega_Li_AM;
      double x_li = lithium_content[i];
      lithium_concentration[i] = x_li * c_li_max / max_lithium_content;
    }
    else if (type[i] == CB_type) {
      li_mols[i] = volume / Omega_Li_CB;
      double x_li = lithium_content[i];
      lithium_concentration[i] = x_li * c_li_max / max_lithium_content;
    }
    else if (type[i] == LM_type) {
      // LM is constant lithium source
      li_mols[i] = volume / Omega_Li_LM;
      lithium_content[i] = max_lithium_content;
      lithium_concentration[i] = c_li_max;
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
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusion::post_force(int vflag)
{
  update_lithium_content();
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusion::updatePtrs()
{
  lithium_content = fix_lithium_content->vector_atom;
  lithium_concentration = fix_lithium_concentration->vector_atom;
  if(fix_current_SE_Li) current_SE_Li = fix_current_SE_Li->vector_atom;
  diffusion_coefficient = fix_diffusion_coefficient->vector_atom;
  lithium_flux = fix_lithium_flux->vector_atom;
  li_mols = fix_li_mols->vector_atom;
}

/* ---------------------------------------------------------------------- */

double FixLithiumDiffusion::get_diffusion_coefficient(int type_i, int type_j)
{
  // Returns the appropriate diffusion coefficient for particle pair
  if (type_i == AM_type) {
    if (type_j == AM_type) return D_AM_AM;
    if (type_j == CB_type) return D_AM_CB;
    if (type_j == LM_type) return D_AM_LM;
  }
  else if (type_i == CB_type) {
    if (type_j == AM_type) return D_CB_AM;
    if (type_j == CB_type) return D_CB_CB;
    if (type_j == LM_type) return D_CB_LM;
  }
  return 0.0;
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
  
  // Reset flux for AM and CB particles only
  for (i = 0; i < nlocal; i++) {
    lithium_flux[i] = 0.0;
  }
  
  // Calculate diffusion flux
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    int type_i = type[i];
    
    // Only process AM and CB particles (LM is constant source)
    if (type_i != AM_type && type_i != CB_type) continue;
    
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      
      int type_j = type[j];
      
      // Check if neighbor is AM, CB, or LM
      if (type_j != AM_type && type_j != CB_type && type_j != LM_type) continue;
      
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      r = sqrt(rsq); // um
      double r_SI = r * 1.0e-6;  // Convert to m
      
      if (r < (radius[i] + radius[j])) {
        double contact_area = calculate_contact_area(i, j);
        
        // Get concentration of neighbor (LM always at c_li_max)
        double c_j = (type_j == LM_type) ? c_li_max : lithium_concentration[j];
        double c_i = lithium_concentration[i];
        
        // Get appropriate diffusion coefficient
        double D_ij = get_diffusion_coefficient(type_i, type_j);
        
        // Diffusion flux: A * D * (c_j - c_i) / |P_j - P_i|
        double c_diff = c_j - c_i;
        double flux = contact_area * D_ij * c_diff / r_SI; // mol/s
        
        lithium_flux[i] += flux;
      }
    }
    
    // Subtract current contribution from SE (if present)
    // -sum(A * i_SE / F)
    if (current_SE_Li) {
      lithium_flux[i] -= current_SE_Li[i] / F;
    }
  }
  
  // Update Notes: The scalling factor needs to be toned down initially to avoid large jumps in lithium content
  // When Lithium content is at least 0.7 then you can use the max scaling factor of 2e16
  // The scaling factor is based on a timestep of 5e-6 us in main, so it simulates 0.1s of EC in 1 run
  double s_factor = 2.0e10; // For 0.1s / 5e-12s
  
  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    int type_i = type[i];
    
    // Skip LM particles - they maintain constant concentration
    if (type_i == LM_type) {
      lithium_concentration[i] = c_li_max;
      lithium_content[i] = max_lithium_content;
      continue;
    }
    
    if (type_i == AM_type || type_i == CB_type) {
      // dn_Li/dt = flux
      lithium_content[i] += (lithium_flux[i] * dt * s_factor) / li_mols[i];
      li_mols[i] += (lithium_flux[i] * dt * s_factor);
 
      // Bounds check
      if (lithium_content[i] < initial_lithium_content) 
        lithium_content[i] = initial_lithium_content;
      if (lithium_content[i] > target_lithium_content) 
        lithium_content[i] = target_lithium_content;
      
      // Update concentration
      lithium_concentration[i] = lithium_content[i] * c_li_max / max_lithium_content;
    }
  }
  
  fix_lithium_content->do_forward_comm();
  fix_lithium_concentration->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

double FixLithiumDiffusion::calculate_contact_area(int i, int j)
{
  double *radius = atom->radius;
  double **x = atom->x;
  double length_conversion = 1.0e-6;
  double delx = (x[i][0] - x[j][0]) * length_conversion;
  double dely = (x[i][1] - x[j][1]) * length_conversion;
  double delz = (x[i][2] - x[j][2]) * length_conversion;
  double r = sqrt(delx*delx + dely*dely + delz*delz);
  double radi = radius[i] * length_conversion;
  double radj = radius[j] * length_conversion;

  double radsum = radi + radj;
  
  if (r >= radsum) return 0.0;
  
  double delta_n = radsum - r;
  double reff = radi * radj / radsum;
  return M_PI * delta_n * reff; // m²
}

/* ---------------------------------------------------------------------- */

double FixLithiumDiffusion::compute_scalar()
{
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double sum = 0.0;
  int count = 0;
  
  for (int i = 0; i < nlocal; i++) {
    if (type[i] == AM_type || type[i] == CB_type) {
      sum += fabs(lithium_flux[i]);
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