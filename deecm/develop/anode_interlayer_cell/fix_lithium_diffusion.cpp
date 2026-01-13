/* ----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine
    Contributing author: Joseph Vazquez Mercado, RIT 2025
    Copyright 2024-     DCS Computing GmbH, Linz
    
    Updated: Concentration + Potential driven diffusion
    - Same-type (AM-AM, CB-CB): concentration-driven with D
    - Cross-type (AM-CB, AM-LM, CB-AM, CB-LM): potential-driven with M
    - M = (c_avg * D * F) / (R * T)
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
  R(8.314),                      // J/(mol·K)
  T(303.0),                      // K (default, can be overridden)
  F(96485.0),                    // C/mol
  c_li_max_LM(77101.002),        // Lithium Metal (mol/m³)
  c_li_max_AM(115917.8),         // Active Material (mol/m³)
  c_li_max_CB(23332.2),          // Carbon Black (mol/m³)
  V_exp_max_AM(3.36),
  V_exp_max_CB(1.13),
  V_exp_max_LM(0.0),
  Omega_Li_LM(12.97e-6),         // m³/mol
  initial_lithium_content(0.0),
  target_lithium_content(1.0),
  max_lithium_content(1.0),
  // Diffusion coefficients for same-type (m²/s)
  D_AM_AM(1.0e-15),
  D_CB_CB(1.0e-15),
  // Diffusion coefficients for cross-type mobility calculation (m²/s)
  D_AM_CB(1.0e-15),
  D_AM_LM(1.0e-15),
  D_CB_AM(1.0e-15),
  D_CB_LM(1.0e-15),
  s_factor(2.0e10),
  lithium_content(NULL),
  lithium_concentration(NULL),
  equilibrium_potential(NULL),
  current_SE_Li(NULL),
  diffusion_coefficient(NULL),
  lithium_flux(NULL),
  li_mols(NULL),
  initial_volume(NULL),
  fix_initial_volume(NULL),
  fix_lithium_content(NULL),
  fix_lithium_concentration(NULL),
  fix_equilibrium_potential(NULL),
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
    } else if (strcmp(arg[iarg],"temperature") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      T = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"c_li_max_LM") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      c_li_max_LM = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"c_li_max_AM") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      c_li_max_AM = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"c_li_max_CB") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      c_li_max_CB = force->numeric(FLERR,arg[iarg+1]);
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
    } else if (strcmp(arg[iarg],"s_factor") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion command");
      s_factor = force->numeric(FLERR,arg[iarg+1]);
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
  fix_diffusion_coefficient = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("diffusionCoefficient","property/atom","scalar",0,0,style,false));
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
  fix_lithium_flux = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("lithiumFlux","property/atom","scalar",0,0,style,false));
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
  fix_li_mols = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("lithiumMols","property/atom","scalar",0,0,style,false));
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
  fix_lithium_concentration = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("lithiumConcentration","property/atom","scalar",0,0,style,false));
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

  // Register property/atom for initial volume
  fix_initial_volume = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("initialVolume","property/atom","scalar",0,0,style,false));
  if(!fix_initial_volume) {
    const char* fixarg[10];
    fixarg[0]="initialVolume";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="initialVolume";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_initial_volume = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
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
  fix_lithium_content = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("lithiumContent","property/atom","scalar",0,0,style));
  if(!fix_lithium_content)
    error->all(FLERR,"Fix lithium_diffusion requires lithiumContent property");
    
  fix_lithium_concentration = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("lithiumConcentration","property/atom","scalar",0,0,style));
  if(!fix_lithium_concentration)
    error->all(FLERR,"Fix lithium_diffusion requires lithiumConcentration property");
  
  fix_initial_volume = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("initialVolume","property/atom","scalar",0,0,style));
  if(!fix_initial_volume)
    error->all(FLERR,"Could not find initialVolume property/atom fix");
  
  // Equilibrium potential is required for potential-driven diffusion
  fix_equilibrium_potential = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("equilibriumPotential","property/atom","scalar",0,0,style));
  if(!fix_equilibrium_potential)
    error->all(FLERR,"Fix lithium_diffusion requires equilibriumPotential property");
  
  fix_current_SE_Li = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("currentSELi","property/atom","scalar",0,0,style,false));
    
  fix_diffusion_coefficient = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("diffusionCoefficient","property/atom","scalar",0,0,style));
  fix_lithium_flux = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("lithiumFlux","property/atom","scalar",0,0,style));
  fix_li_mols = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("lithiumMols","property/atom","scalar",0,0,style));
  
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
    double c_max_i = get_c_li_max(type[i]);
    
    initial_volume[i] = volume;

    if (type[i] == AM_type) {
      double x_li = lithium_content[i];
      li_mols[i] = x_li * c_max_i * volume / max_lithium_content;
      lithium_concentration[i] = x_li * c_max_i / max_lithium_content;
    }
    else if (type[i] == CB_type) {
      double x_li = lithium_content[i];
      li_mols[i] = x_li * c_max_i * volume / max_lithium_content;
      lithium_concentration[i] = x_li * c_max_i / max_lithium_content;
    }
    else if (type[i] == LM_type) {
      li_mols[i] = volume / Omega_Li_LM;
      lithium_content[i] = max_lithium_content;
      lithium_concentration[i] = c_max_i;
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
  // CRITICAL: Ensure ghost atoms have current values BEFORE calculating fluxes
  fix_lithium_content->do_forward_comm();
  fix_lithium_concentration->do_forward_comm();
  fix_equilibrium_potential->do_forward_comm();

  update_lithium_content();
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusion::updatePtrs()
{
  lithium_content = fix_lithium_content->vector_atom;
  lithium_concentration = fix_lithium_concentration->vector_atom;
  equilibrium_potential = fix_equilibrium_potential->vector_atom;
  if(fix_current_SE_Li) current_SE_Li = fix_current_SE_Li->vector_atom;
  diffusion_coefficient = fix_diffusion_coefficient->vector_atom;
  lithium_flux = fix_lithium_flux->vector_atom;
  li_mols = fix_li_mols->vector_atom;
  initial_volume = fix_initial_volume->vector_atom;
}

/* ---------------------------------------------------------------------- */

double FixLithiumDiffusion::get_diffusion_coefficient(int type_i, int type_j)
{
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

double FixLithiumDiffusion::get_c_li_max(int type)
{
  if (type == AM_type) return c_li_max_AM;
  if (type == CB_type) return c_li_max_CB;
  if (type == LM_type) return c_li_max_LM;
  return 0.0;
}

/* ---------------------------------------------------------------------- */

double FixLithiumDiffusion::get_V_exp_max(int type)
{
  if (type == AM_type) return V_exp_max_AM;
  if (type == CB_type) return V_exp_max_CB;
  if (type == LM_type) return V_exp_max_LM;
  return 0.0;
}

/* ---------------------------------------------------------------------- */

double FixLithiumDiffusion::compute_mobility(double c_i, double c_j, double D_ij)
{
  // M = (c_avg * D * F) / (R * T)
  double c_avg = 0.5 * (c_i + c_j);
  return (c_avg * D_ij * F) / (R * T);
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
  
  // Reset flux for AM and CB particles
  for (i = 0; i < nlocal; i++) {
    lithium_flux[i] = 0.0;
  }
  
  // Calculate diffusion flux based on new equations
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    int type_i = type[i];
    
    // Only process AM and CB particles (LM is constant source)
    if (type_i != AM_type && type_i != CB_type) continue;
    
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    
    double c_i = lithium_concentration[i];
    double U_i = equilibrium_potential[i];
    
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
        
        // Get neighbor properties
        double c_j = (type_j == LM_type) ? c_li_max_LM : lithium_concentration[j];
        double U_j = equilibrium_potential[j];
        
        double D_ij = get_diffusion_coefficient(type_i, type_j);
        double flux = 0.0;
        
        // Determine flux type based on particle types
        if (type_i == AM_type) {
          if (type_j == AM_type) {
            // AM-AM: concentration-driven
            // A * D_AM-AM * (c_j - c_i) / |P_j - P_i|
            flux = contact_area * D_ij * (c_j - c_i) / r_SI;
          } else {
            // AM-CB or AM-LM: potential-driven with mobility
            // A * M * (U_i - U_j) / |P_j - P_i|
            double M_ij = compute_mobility(c_i, c_j, D_ij);
            flux = contact_area * M_ij * (U_i - U_j) / r_SI;
          }
        } else if (type_i == CB_type) {
          if (type_j == CB_type) {
            // CB-CB: concentration-driven
            // A * D_CB-CB * (c_j - c_i) / |P_j - P_i|
            flux = contact_area * D_ij * (c_j - c_i) / r_SI;
          } else {
            // CB-AM or CB-LM: potential-driven with mobility
            // A * M * (U_i - U_j) / |P_j - P_i|
            double M_ij = compute_mobility(c_i, c_j, D_ij);
            flux = contact_area * M_ij * (U_i - U_j) / r_SI;
          }
        }
        
        lithium_flux[i] += flux;
      }
    }
    
    // Subtract current contribution from SE (if present)
    // -sum(A * i_SE / F)
    if (current_SE_Li) {
      lithium_flux[i] -= current_SE_Li[i] / F;
    }
  }
  
  // Scaling factor for time acceleration
  // double s_factor = 2.0e10;
  
  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    int type_i = type[i];
    
    // Skip LM particles - they maintain constant concentration
    if (type_i == LM_type) {
      lithium_concentration[i] = c_li_max_LM;
      lithium_content[i] = max_lithium_content;
      continue;
    }
    
    if (type_i == AM_type || type_i == CB_type) {
      double c_max_current = get_c_li_max(type_i);
      double V_exp_max_current = get_V_exp_max(type_i);
      
      // dn_Li/dt = flux -> update lithium content
      lithium_content[i] += (lithium_flux[i] * dt * s_factor) * max_lithium_content / 
                           (c_max_current * initial_volume[i] * V_exp_max_current);
 
      // Bounds check
      if (lithium_content[i] < initial_lithium_content)
        lithium_content[i] = initial_lithium_content;
      if (lithium_content[i] > target_lithium_content)
        lithium_content[i] = target_lithium_content;
      
      // Update Number of Moles in Particle
      li_mols[i] = lithium_content[i] * c_max_current * initial_volume[i] * 
                   V_exp_max_current / max_lithium_content;
      
      // Update concentration
      lithium_concentration[i] = lithium_content[i] * c_max_current / max_lithium_content;
    }
  }
  
  fix_lithium_content->do_forward_comm();
  fix_lithium_concentration->do_forward_comm();
  fix_li_mols->do_forward_comm();  // ADD THIS - needed by radius_expansion
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