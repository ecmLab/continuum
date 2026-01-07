/* ----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine
    http://www.dcs-computing.com, office@dcs-computing.com
    LIGGGHTS® is part of CFDEM®project: http://www.liggghts.com | http://www.cfdem.com

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    Created: Joseph Vazquez Mercado, RIT 2025
    Copyright 2024-     DCS Computing GmbH, Linz
    Notes: Metal-Carbon system lithium diffusion
           - Carbon particles act as fully lithiated pathways
           - Metal particles receive lithium via concentration-based diffusion
           - Supports infinite lithium source region (e.g., z=0 plane)
------------------------------------------------------------------------- */

#include "fix_metal_carbon_diffusion.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "error.h"
#include "force.h"
#include "fix_property_atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "domain.h"
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMetalCarbonDiffusion::FixMetalCarbonDiffusion(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  c_li_max(76900.0),           // mol/m³ max Li concentration in fully lithiated carbon
  c_li_source(76900.0),        // mol/m³ Li concentration at source region
  D_Li_C(1.0e-14),             // m²/s Li diffusion in carbon
  D_Li_metal(1.0e-14),         // m²/s Li diffusion in metal
  source_type(SOURCE_NONE),
  source_coord(0.0),
  source_thickness(0.01),      // μm thickness of source region
  lithium_content(NULL),
  lithium_concentration(NULL),
  diffusion_coefficient(NULL),
  lithium_flux(NULL),
  li_mols(NULL),
  metal_mols(NULL),
  fix_lithium_content(NULL),
  fix_lithium_concentration(NULL),
  fix_diffusion_coefficient(NULL),
  fix_lithium_flux(NULL),
  fix_li_mols(NULL),
  fix_metal_mols(NULL),
  metal_type(1),
  carbon_type(2),
  list(NULL)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix metal_carbon_diffusion command");

  // Parse arguments
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"metal_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix metal_carbon_diffusion command");
      metal_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"carbon_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix metal_carbon_diffusion command");
      carbon_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"c_li_max") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix metal_carbon_diffusion command");
      c_li_max = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"D_Li_C") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix metal_carbon_diffusion command");
      D_Li_C = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"D_Li_metal") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix metal_carbon_diffusion command");
      D_Li_metal = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"source") == 0) {
      // source zlo/zhi/xlo/xhi/ylo/yhi coord thickness c_source
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix metal_carbon_diffusion command: source requires 4 args");
      if (strcmp(arg[iarg+1],"zlo") == 0) source_type = SOURCE_ZLO;
      else if (strcmp(arg[iarg+1],"zhi") == 0) source_type = SOURCE_ZHI;
      else if (strcmp(arg[iarg+1],"xlo") == 0) source_type = SOURCE_XLO;
      else if (strcmp(arg[iarg+1],"xhi") == 0) source_type = SOURCE_XHI;
      else if (strcmp(arg[iarg+1],"ylo") == 0) source_type = SOURCE_YLO;
      else if (strcmp(arg[iarg+1],"yhi") == 0) source_type = SOURCE_YHI;
      else error->all(FLERR,"Unknown source type for fix metal_carbon_diffusion");
      source_coord = force->numeric(FLERR,arg[iarg+2]);
      source_thickness = force->numeric(FLERR,arg[iarg+3]);
      c_li_source = force->numeric(FLERR,arg[iarg+4]);
      iarg += 5;
    } else error->all(FLERR,"Illegal fix metal_carbon_diffusion command");
  }

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;
}

/* ---------------------------------------------------------------------- */

void FixMetalCarbonDiffusion::post_create()
{
  // Register property/atom for diffusion coefficient
  fix_diffusion_coefficient = static_cast<FixPropertyAtom*>(modify->find_fix_property("diffusionCoefficient","property/atom","scalar",0,0,style,false));
  if(!fix_diffusion_coefficient) {
    const char* fixarg[9];
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
    const char* fixarg[9];
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
    const char* fixarg[9];
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

  // Register property/atom for metal mols
  fix_metal_mols = static_cast<FixPropertyAtom*>(modify->find_fix_property("metalMols","property/atom","scalar",0,0,style,false));
  if(!fix_metal_mols) {
    const char* fixarg[9];
    fixarg[0]="metalMols";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="metalMols";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_metal_mols = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  // Register property/atom for lithium concentration
  fix_lithium_concentration = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumConcentration","property/atom","scalar",0,0,style,false));
  if(!fix_lithium_concentration) {
    const char* fixarg[9];
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

  // Register property/atom for lithium content (Li/Metal ratio)
  fix_lithium_content = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumContent","property/atom","scalar",0,0,style,false));
  if(!fix_lithium_content) {
    const char* fixarg[9];
    fixarg[0]="lithiumContent";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="lithiumContent";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_lithium_content = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }
}

/* ---------------------------------------------------------------------- */

FixMetalCarbonDiffusion::~FixMetalCarbonDiffusion()
{
}

/* ---------------------------------------------------------------------- */

int FixMetalCarbonDiffusion::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMetalCarbonDiffusion::init()
{
  // Find required fixes
  fix_lithium_content = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumContent","property/atom","scalar",0,0,style));
  if(!fix_lithium_content)
    error->all(FLERR,"Fix metal_carbon_diffusion requires lithiumContent property");
    
  fix_lithium_concentration = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumConcentration","property/atom","scalar",0,0,style));
  if(!fix_lithium_concentration)
    error->all(FLERR,"Fix metal_carbon_diffusion requires lithiumConcentration property");
    
  fix_diffusion_coefficient = static_cast<FixPropertyAtom*>(modify->find_fix_property("diffusionCoefficient","property/atom","scalar",0,0,style));
  fix_lithium_flux = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumFlux","property/atom","scalar",0,0,style));
  fix_li_mols = static_cast<FixPropertyAtom*>(modify->find_fix_property("lithiumMols","property/atom","scalar",0,0,style));
  fix_metal_mols = static_cast<FixPropertyAtom*>(modify->find_fix_property("metalMols","property/atom","scalar",0,0,style));
  
  if(!fix_diffusion_coefficient || !fix_lithium_flux || !fix_li_mols || !fix_metal_mols)
    error->all(FLERR,"Could not find required property/atom fixes");

  // Request neighbor list
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixMetalCarbonDiffusion::setup(int vflag)
{
  updatePtrs();
  
  double *radius = atom->radius;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  int *mask = atom->mask;
  
  double Omega_metal = 10.96e-6;  // m³/mol molar volume of metal (e.g., Ag)
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double radius_m = radius[i] * 1.0e-6;  // Convert μm to m
      double V_particle = (4.0/3.0) * M_PI * radius_m * radius_m * radius_m;  // m³
      
      if (type[i] == metal_type) {
        // Metal particles: initially no lithium, calculate metal mols from volume
        metal_mols[i] = V_particle / Omega_metal;  // mols of metal
        li_mols[i] = 0.0;  // Initially no lithium
        lithium_content[i] = 0.0;  // Li/Metal ratio = 0
        lithium_concentration[i] = 0.0;  // mol/m³
        diffusion_coefficient[i] = D_Li_metal;
        
      } else if (type[i] == carbon_type) {
        // Carbon particles: always fully lithiated (act as Li pathway)
        metal_mols[i] = 0.0;  // No metal in carbon
        lithium_concentration[i] = c_li_max;  // Fully lithiated
        lithium_content[i] = 1.0;  // Fully lithiated ratio
        li_mols[i] = c_li_max * V_particle;  // mols of Li in carbon
        diffusion_coefficient[i] = D_Li_C;
      }
      
      // Check if particle is in source region and set concentration
      if (is_in_source_region(x[i])) {
        lithium_concentration[i] = c_li_source;
        lithium_content[i] = 1.0;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMetalCarbonDiffusion::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixMetalCarbonDiffusion::pre_force(int vflag)
{
  updatePtrs();
  
  // Maintain carbon particles as fully lithiated (they are pathways)
  int nlocal = atom->nlocal;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      // Carbon always fully lithiated
      if (type[i] == carbon_type) {
        lithium_concentration[i] = c_li_max;
        lithium_content[i] = 1.0;
      }
      // Particles in source region maintain source concentration
      if (is_in_source_region(x[i])) {
        lithium_concentration[i] = c_li_source;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMetalCarbonDiffusion::post_force(int vflag)
{
  update_lithium_content();
}

/* ---------------------------------------------------------------------- */

void FixMetalCarbonDiffusion::updatePtrs()
{
  lithium_content = fix_lithium_content->vector_atom;
  lithium_concentration = fix_lithium_concentration->vector_atom;
  diffusion_coefficient = fix_diffusion_coefficient->vector_atom;
  lithium_flux = fix_lithium_flux->vector_atom;
  li_mols = fix_li_mols->vector_atom;
  metal_mols = fix_metal_mols->vector_atom;
}

/* ---------------------------------------------------------------------- */

bool FixMetalCarbonDiffusion::is_in_source_region(double *pos)
{
  if (source_type == SOURCE_NONE) return false;
  
  double coord;
  switch (source_type) {
    case SOURCE_ZLO:
      coord = pos[2];
      return (coord >= source_coord && coord <= source_coord + source_thickness);
    case SOURCE_ZHI:
      coord = pos[2];
      return (coord <= source_coord && coord >= source_coord - source_thickness);
    case SOURCE_XLO:
      coord = pos[0];
      return (coord >= source_coord && coord <= source_coord + source_thickness);
    case SOURCE_XHI:
      coord = pos[0];
      return (coord <= source_coord && coord >= source_coord - source_thickness);
    case SOURCE_YLO:
      coord = pos[1];
      return (coord >= source_coord && coord <= source_coord + source_thickness);
    case SOURCE_YHI:
      coord = pos[1];
      return (coord <= source_coord && coord >= source_coord - source_thickness);
    default:
      return false;
  }
}

/* ---------------------------------------------------------------------- */

void FixMetalCarbonDiffusion::update_lithium_content()
{
  int i,j,ii,jj,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,r;
  
  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt = update->dt * 1.0e-6;  // Convert μs to s

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // Reset flux for metal particles only (carbon maintains constant concentration)
  for (i = 0; i < nlocal; i++) {
    if (type[i] == metal_type) {
      lithium_flux[i] = 0.0;
    }
  }
  
  // Calculate diffusion flux between particles
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    
    // Only calculate flux INTO metal particles
    if (type[i] != metal_type) continue;
    
    // Skip if particle is in source region (infinite source)
    if (is_in_source_region(x[i])) continue;
    
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      
      // Diffusion from carbon to metal, or metal to metal
      if (type[j] == carbon_type || type[j] == metal_type) {
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        r = sqrt(rsq);  // μm
        double r_SI = r * 1.0e-6;  // Convert to m
        
        if (r < (radius[i] + radius[j])) {
          // Calculate contact area
          double contact_area = calculate_contact_area(i, j);
          
          // Concentration gradient (from j to i)
          double c_diff = lithium_concentration[j] - lithium_concentration[i];  // mol/m³
          
          // Use harmonic mean of diffusion coefficients at interface
          double D_eff;
          if (type[j] == carbon_type) {
            // Carbon-metal interface
            D_eff = 2.0 * D_Li_C * D_Li_metal / (D_Li_C + D_Li_metal);
          } else {
            // Metal-metal interface
            D_eff = D_Li_metal;
          }
          
          // Fick's first law: flux = -D * A * (dc/dx)
          // Positive flux means Li flowing into particle i
          double flux = contact_area * D_eff * c_diff / r_SI;  // mol/s
          
          lithium_flux[i] += flux;
        }
      }
    }
    
    // Add flux from source region if particle is adjacent to it
    if (source_type != SOURCE_NONE) {
      double dist_to_source = distance_to_source(x[i]);
      if (dist_to_source < radius[i] * 1.0e-6) {  // Within contact distance
        // Approximate contact area with source plane
        double radius_m = radius[i] * 1.0e-6;
        double overlap = radius_m - dist_to_source;
        if (overlap > 0) {
          double source_contact_area = M_PI * overlap * radius_m;
          double c_diff = c_li_source - lithium_concentration[i];
          double flux = source_contact_area * D_Li_metal * c_diff / radius_m;
          lithium_flux[i] += flux;
        }
      }
    }
  }
  
  // Update lithium content for metal particles based on flux
  // Scaling factor to accelerate diffusion (adjust based on timestep)
  double s_factor = 2.0e10;  // Scale factor for time acceleration
  
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && type[i] == metal_type) {
      // Skip source region particles
      if (is_in_source_region(x[i])) continue;
      
      // Update lithium mols: dn_Li/dt = flux
      double delta_li = lithium_flux[i] * dt * s_factor;
      li_mols[i] += delta_li;
      
      // Ensure non-negative lithium content
      if (li_mols[i] < 0.0) li_mols[i] = 0.0;
      
      // Update lithium content ratio (Li/Metal)
      if (metal_mols[i] > 0.0) {
        lithium_content[i] = li_mols[i] / metal_mols[i];
      }
      
      // Update concentration based on current particle volume
      double radius_m = radius[i] * 1.0e-6;
      double V_particle = (4.0/3.0) * M_PI * radius_m * radius_m * radius_m;
      if (V_particle > 0.0) {
        lithium_concentration[i] = li_mols[i] / V_particle;
      }
    }
  }
  
  // Forward communication to update ghost particles
  fix_lithium_content->do_forward_comm();
  fix_lithium_concentration->do_forward_comm();
  fix_li_mols->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

double FixMetalCarbonDiffusion::distance_to_source(double *pos)
{
  if (source_type == SOURCE_NONE) return 1.0e30;
  
  switch (source_type) {
    case SOURCE_ZLO:
    case SOURCE_ZHI:
      return fabs(pos[2] - source_coord);
    case SOURCE_XLO:
    case SOURCE_XHI:
      return fabs(pos[0] - source_coord);
    case SOURCE_YLO:
    case SOURCE_YHI:
      return fabs(pos[1] - source_coord);
    default:
      return 1.0e30;
  }
}

/* ---------------------------------------------------------------------- */

double FixMetalCarbonDiffusion::calculate_contact_area(int i, int j)
{
  double *radius = atom->radius;
  double **x = atom->x;
  double length_conversion = 1.0e-6;  // μm to m
  
  double delx = (x[i][0] - x[j][0]) * length_conversion;
  double dely = (x[i][1] - x[j][1]) * length_conversion;
  double delz = (x[i][2] - x[j][2]) * length_conversion;
  double r = sqrt(delx*delx + dely*dely + delz*delz);  // m
  
  double radi = radius[i] * length_conversion;  // m
  double radj = radius[j] * length_conversion;  // m
  double radsum = radi + radj;
  
  if (r >= radsum) return 0.0;
  
  // Hertzian contact area approximation
  double delta_n = radsum - r;  // Overlap
  double reff = radi * radj / radsum;  // Effective radius
  return M_PI * delta_n * reff;  // m²
}

/* ---------------------------------------------------------------------- */

double FixMetalCarbonDiffusion::compute_scalar()
{
  // Return average lithium flux magnitude
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double sum = 0.0;
  int count = 0;
  
  for (int i = 0; i < nlocal; i++) {
    if (type[i] == metal_type) {
      sum += fabs(lithium_flux[i]);
      count++;
    }
  }
  
  double all_sum;
  int all_count;
  MPI_Allreduce(&sum,&all_sum,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&count,&all_count,1,MPI_INT,MPI_SUM,world);
  
  if (all_count > 0) return all_sum / all_count;
  return 0.0;
}