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
    Created by: Joseph Vazquez
    Copyright 2024-     DCS Computing GmbH, Linz

    Notes: 
    - FOR CONTROLED VOLTAGE
    - Solves both electronic and electrolyte potentials
------------------------------------------------------------------------- */

#include "fix_battery_eis.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "error.h"
#include "force.h"
#include "fix_property_atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "group.h"
#include "domain.h"
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL 1e-15
#define MAX_EXP_ARG 20.0  // Prevent exp() overflow

/* ---------------------------------------------------------------------- */

void FixBatteryEIS::post_create()
{
  // Register property/atom for electrolyte potential
  fix_phi_el = static_cast<FixPropertyAtom*>(modify->find_fix_property("electrolytePotential","property/atom","scalar",0,0,style,false));
  if(!fix_phi_el) {
    const char* fixarg[10];
    fixarg[0]="electrolytePotential";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="electrolytePotential";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_phi_el = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  // Register property/atom for old electrolyte potential (convergence check)
  fix_phi_el_old = static_cast<FixPropertyAtom*>(modify->find_fix_property("electrolytePotentialOld","property/atom","scalar",0,0,style,false));
  if(!fix_phi_el_old) {
    const char* fixarg[10];
    fixarg[0]="electrolytePotentialOld";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="electrolytePotentialOld";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_phi_el_old = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  // Register property/atom for electronic potential
  fix_phi_ed = static_cast<FixPropertyAtom*>(modify->find_fix_property("electronicPotential","property/atom","scalar",0,0,style,false));
  if(!fix_phi_ed) {
    const char* fixarg[10];
    fixarg[0]="electronicPotential";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="electronicPotential";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_phi_ed = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  // Register property/atom for old electronic potential (convergence check)
  fix_phi_ed_old = static_cast<FixPropertyAtom*>(modify->find_fix_property("electronicPotentialOld","property/atom","scalar",0,0,style,false));
  if(!fix_phi_ed_old) {
    const char* fixarg[10];
    fixarg[0]="electronicPotentialOld";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="electronicPotentialOld";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_phi_ed_old = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  // Register property/atom for current density Li-SE
  fix_current_Li_SE = static_cast<FixPropertyAtom*>(modify->find_fix_property("currentLiSE","property/atom","scalar",0,0,style,false));
  if(!fix_current_Li_SE) {
    const char* fixarg[10];
    fixarg[0]="currentLiSE";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="currentLiSE";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_current_Li_SE = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  // Register property/atom for hydrostatic stress
  fix_hydrostatic_stress = static_cast<FixPropertyAtom*>(modify->find_fix_property("hydrostaticStress","property/atom","scalar",0,0,style,false));
  if(!fix_hydrostatic_stress) {
    const char* fixarg[10];
    fixarg[0]="hydrostaticStress";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="hydrostaticStress";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_hydrostatic_stress = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  // Register property/atom for initialization flag
  fix_init_flag = static_cast<FixPropertyAtom*>(modify->find_fix_property("batteryInitFlag","property/atom","scalar",0,0,style,false));
  if(!fix_init_flag) {
    const char* fixarg[10];
    fixarg[0]="batteryInitFlag";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="batteryInitFlag";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";  // 0 = not initialized
    fix_init_flag = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }
}

/* ---------------------------------------------------------------------- */

FixBatteryEIS::FixBatteryEIS(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  omega(1.9),
  tolerance(1e-6),
  max_iterations(10000),
  current_iteration(0),
  convergence_error(0.0),
  R(8.31), // J*mol^(-1)*K^(-1)
  T(303.0), // K
  F(96485.0), // C/mol
  sigma_el(0.05), // Electrolyte conductivity (S/m) - for SE
  sigma_ed_SE(200.0), // Electronic conductivity for CBD/SE (S/m)
  sigma_ed_CC(300.0), // Electronic conductivity for CC (S/m)
  alpha_a(0.5), //unitless
  alpha_c(0.5),
  BC_bottom_type(2),      // Default: CC1 Refrence
  BC_top_type(3),         // Default: CC2 Anode
  phi_el_BC_bottom(0.0),  // Default: 0V Initially for electrolyte
  phi_el_BC_top(0.01),     // Default: 10mV Always on Anode for electrolyte
  phi_ed_BC_anode(0.0),   // Electronic potential at anode = 0V
  phi_el(NULL),
  phi_el_old(NULL),
  phi_ed(NULL),
  phi_ed_old(NULL),
  current_Li_SE(NULL),
  hydrostatic_stress(NULL),
  fix_phi_el(NULL),
  fix_phi_el_old(NULL),
  fix_phi_ed(NULL),
  fix_phi_ed_old(NULL),
  fix_current_Li_SE(NULL),
  fix_hydrostatic_stress(NULL),
  fix_init_flag(NULL),
  groupbit_SE(0),
  SE_type(1),
  list(NULL)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix battery/eis command");

  // Parse arguments
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"omega") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix battery/eis command");
      omega = force->numeric(FLERR,arg[iarg+1]);
      if (omega <= 0.0 || omega >= 2.0)
        error->all(FLERR,"EIS omega must be between 0 and 2");
      iarg += 2;
    } else if (strcmp(arg[iarg],"tolerance") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix battery/eis command");
      tolerance = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"max_iter") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix battery/eis command");
      max_iterations = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"temperature") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix battery/eis command");
      T = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"BC_types") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix battery/eis command");
      BC_bottom_type = force->inumeric(FLERR,arg[iarg+1]);
      BC_top_type = force->inumeric(FLERR,arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"BC_potentials") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix battery/eis command");
      phi_el_BC_bottom = force->numeric(FLERR,arg[iarg+1]);
      phi_el_BC_top = force->numeric(FLERR,arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"conductivity") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix battery/eis command");
      sigma_el = force->numeric(FLERR,arg[iarg+1]);      // Electrolyte
      iarg += 2;
    } else if (strcmp(arg[iarg],"SE_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix battery/eis command");
      SE_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix battery/eis command");
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

FixBatteryEIS::~FixBatteryEIS()
{
}

/* ---------------------------------------------------------------------- */

int FixBatteryEIS::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBatteryEIS::init()
{
  // Check newton pair
  if (force->newton_pair == 1)
    error->all(FLERR,"Fix battery/eis requires newton pair off");

  // Get all property fixes
  fix_phi_el = static_cast<FixPropertyAtom*>(modify->find_fix_property("electrolytePotential","property/atom","scalar",0,0,style));
  fix_phi_el_old = static_cast<FixPropertyAtom*>(modify->find_fix_property("electrolytePotentialOld","property/atom","scalar",0,0,style));
  fix_phi_ed = static_cast<FixPropertyAtom*>(modify->find_fix_property("electronicPotential","property/atom","scalar",0,0,style));
  fix_phi_ed_old = static_cast<FixPropertyAtom*>(modify->find_fix_property("electronicPotentialOld","property/atom","scalar",0,0,style));
  fix_current_Li_SE = static_cast<FixPropertyAtom*>(modify->find_fix_property("currentLiSE","property/atom","scalar",0,0,style));
  fix_hydrostatic_stress = static_cast<FixPropertyAtom*>(modify->find_fix_property("hydrostaticStress","property/atom","scalar",0,0,style));
  fix_init_flag = static_cast<FixPropertyAtom*>(modify->find_fix_property("batteryInitFlag","property/atom","scalar",0,0,style));

  if(!fix_phi_el || !fix_phi_el_old || !fix_phi_ed || !fix_phi_ed_old || !fix_current_Li_SE || !fix_hydrostatic_stress || !fix_init_flag)
    error->all(FLERR,"Could not find required property/atom fixes");

  // Validate particle types
  if (SE_type < 1 || SE_type > atom->ntypes)
    error->all(FLERR,"Invalid particle types for battery/eis");
    
  if (BC_bottom_type < 1 || BC_bottom_type > atom->ntypes || BC_top_type < 1 || BC_top_type > atom->ntypes)
    error->all(FLERR,"Invalid boundary condition types for battery/eis");

  // Request neighbor list
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixBatteryEIS::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixBatteryEIS::setup(int vflag)
{
  updatePtrs();
  
  int nlocal = atom->nlocal;
  int *type = atom->type;
  int *mask = atom->mask;
  
  double *init_flag = fix_init_flag->vector_atom;

  // Check if potentials have already been initialized
  bool already_initialized = false;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && init_flag[i] > 0.5) {
      already_initialized = true;
      break;
    }
  }
  
  // Only initialize if not already done
  if (!already_initialized) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        // Initialize both potentials based on particle type
        if (type[i] == SE_type) {
          phi_el[i] = 0.0;
          phi_el_old[i] = 0.0;
          phi_ed[i] = 0.0;
          phi_ed_old[i] = 0.0;
        } else if (type[i] == BC_bottom_type) {
          phi_el[i] = phi_el_BC_bottom;
          phi_el_old[i] = phi_el_BC_bottom;
          phi_ed[i] = 0.0;  
          phi_ed_old[i] = 0.0;
        } else if (type[i] == BC_top_type) {
          phi_el[i] = phi_el_BC_top;
          phi_el_old[i] = phi_el_BC_top;
          phi_ed[i] = phi_ed_BC_anode;
          phi_ed_old[i] = phi_ed_BC_anode;
        }
        init_flag[i] = 1.0;
      }
    }
  } else {
    // If already initialized, just update the old values
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        phi_el_old[i] = phi_el[i];
        phi_ed_old[i] = phi_ed[i];
      }
    }
  }
  
  // Apply boundary conditions
  apply_boundary_conditions();
}

/* ---------------------------------------------------------------------- */

void FixBatteryEIS::pre_force(int vflag)
{
  updatePtrs();
  
  // Calculate hydrostatic stress for Li particles
  calculate_hydrostatic_stress();
}

/* ---------------------------------------------------------------------- */

void FixBatteryEIS::post_force(int vflag)
{
  // Phase 1: Solve with Li-SE interface currents
  
  current_iteration = 0;
  convergence_error = 1.0;
  
  while (current_iteration < max_iterations && convergence_error > tolerance) {
    solve_eis_iteration();
    convergence_error = check_convergence();
    current_iteration++;
  }
  
  if (comm->me == 0) {
    if (current_iteration >= max_iterations) {
      error->warning(FLERR,"Battery EIS Phase 1 did not converge");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixBatteryEIS::updatePtrs()
{
  phi_el = fix_phi_el->vector_atom;
  phi_el_old = fix_phi_el_old->vector_atom;
  phi_ed = fix_phi_ed->vector_atom;
  phi_ed_old = fix_phi_ed_old->vector_atom;
  current_Li_SE = fix_current_Li_SE->vector_atom;
  hydrostatic_stress = fix_hydrostatic_stress->vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixBatteryEIS::solve_eis_iteration()
{
  int i,j,ii,jj,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,r,r_SI;
  double phi_el_sum,coeff_el_sum,cur_sum,cur_ed;
  double phi_ed_sum,coeff_ed_sum;
  
  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // Store old values for all particles
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      phi_el_old[i] = phi_el[i];
      phi_ed_old[i] = phi_ed[i];
    }
  }

  // Reset current_Li_SE for all particles
  for (int i = 0; i < nlocal; i++) {
    current_Li_SE[i] = 0.0;
  }

  // EIS iteration for both potentials
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    
    // Skip top and bottom boundary particles as they have fixed potentials
    if (type[i] == BC_top_type || type[i] == BC_bottom_type) continue;
    
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    
    phi_el_sum = 0.0;
    coeff_el_sum = 0.0;
    cur_sum = 0.0;
    cur_ed = 0.0;
    phi_ed_sum = 0.0;
    coeff_ed_sum = 0.0;
    
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    // Determine electronic conductivity for particle i
    double sigma_ed_i = sigma_ed_SE;
    // if (type[i] == SE_type) sigma_ed_i = sigma_ed_SE;
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      r = sqrt(rsq);  // micrometers
      r_SI = r * 1.0e-6;  // Convert to m
      
      if (r < (radius[i] + radius[j])) {
        double contact_area = calculate_contact_area(i, j); // m2
        
        if (contact_area > 0.0) {
          // === Electrolyte potential update (SE and Bottom CC particles only) ===
          if (type[i] == SE_type) {
            if (type[j] == SE_type) {
              double conductance = sigma_el * contact_area / r_SI;
              phi_el_sum += conductance * phi_el[j];
              coeff_el_sum += conductance;
            
            } else if (type[j] == BC_top_type || type[j] == BC_bottom_type) {
              double conductance = sigma_el * contact_area / r_SI;
              phi_el_sum += conductance * phi_el[j];
              coeff_el_sum += conductance;
              double i_pq = calculate_current_Li_SE(j, i, phi_el[j], phi_el[i], hydrostatic_stress[j]);
              current_Li_SE[i] += -1 * i_pq * contact_area;
                
              // Add current contribution to SE particle
              cur_sum += i_pq * contact_area;
            }
          }
        }
      }
    }
    
    // Update electrolyte potential for SE particles only
    if ((type[i] == SE_type) && coeff_el_sum > SMALL) {
      double phi_el_new = (phi_el_sum + cur_sum) / coeff_el_sum;
      
      if (!std::isnan(phi_el_new) && !std::isinf(phi_el_new)) {
        phi_el[i] = phi_el_old[i] + omega * (phi_el_new - phi_el_old[i]);
      } else {
        phi_el[i] = phi_el_old[i];
      }
    }
  }
  
  // Apply boundary conditions after each iteration
  apply_boundary_conditions();
  
  // Forward communication for both potentials
  fix_phi_el->do_forward_comm();
  fix_phi_ed->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

void FixBatteryEIS::apply_boundary_conditions()
{
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  
  // Get current simulation time in microseconds (since units are "micro")
  double t = update->dt * update->ntimestep;
  t = t * 1.0e-6; // Convert to seconds

  // Calculate sinusoidal boundary condition: 0.005*sin(t*5e10) + 0.005
  double phi_el_BC_top_sinusoidal = 0.005 * sin(t * 5.0e10) + 0.005;
  // double phi_el_BC_top_sinusoidal = 0.01; // Keeping constant for debugging

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      // Top boundary condition particles (Anode) - fixed potentials
      if (type[i] == BC_top_type) {
        phi_el[i] = phi_el_BC_top_sinusoidal;      // Electrolyte potential not fixed at anode
        phi_el_old[i] = phi_el_BC_top_sinusoidal;
        // phi_ed[i] = phi_ed_BC_anode;    // Fixed electronic potential (0V) Ignoring electric for now
        // phi_ed_old[i] = phi_ed_BC_anode;
      }
      // Bottom boundary condition particles (CC1) - fixed potentials
      else if (type[i] == BC_bottom_type) {
        phi_el[i] = phi_el_BC_bottom;  // Fixed electrolyte potential (0V)
        phi_el_old[i] = phi_el_BC_bottom;
        // phi_ed[i] = 0.0;               // Electronic potential not fixed at CC1
        // phi_ed_old[i] = 0.0;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixBatteryEIS::calculate_hydrostatic_stress()
{
  double **f = atom->f;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && type[i] == BC_bottom_type) {
      // Calculate magnitude of force
      double force_conversion = 1.0e-9;  // LAMMPS force micro nN to N SI
      double fmag = sqrt(f[i][0]*f[i][0] + f[i][1]*f[i][1] + f[i][2]*f[i][2]) * force_conversion;  // N
      
      // Calculate surface area
      double radius_SI = radius[i] * 1.0e-6;  // Convert to m
      double surface_area = 4.0 * M_PI * radius_SI * radius_SI;  // m2
      
      // Hydrostatic stress = force / area
      if (surface_area > 0.0) {
        // hydrostatic_stress[i] = fmag / surface_area;
        hydrostatic_stress[i] = 0.0; // Setting to zero for debugging
      } else {
        hydrostatic_stress[i] = 0.0;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixBatteryEIS::calculate_contact_area(int i, int j)
{
  double *radius = atom->radius;
  double **x = atom->x;
  
  double length_conversion = 1.0e-6;
  double delx = (x[i][0] - x[j][0]) * length_conversion;
  double dely = (x[i][1] - x[j][1]) * length_conversion;
  double delz = (x[i][2] - x[j][2]) * length_conversion;
  double r = sqrt(delx*delx + dely*dely + delz*delz);  // m
  
  double radi = radius[i] * length_conversion;  // m
  double radj = radius[j] * length_conversion;  // m
  double radsum = radi + radj;
  
  if (r >= radsum) return 0.0;
  
  // Overlap-based contact area
  if (r < fmax(radi, radj)) {
    // One sphere inside the other
    double rmin = fmin(radi, radj);
    return M_PI * rmin * rmin;
  } else {
    // Hertz contact area approximation
    double delta_n = radsum - r;
    double reff = radi * radj / radsum;
    // Limit overlap to prevent unrealistic areas
    if (delta_n > 0.1 * fmin(radi, radj)) {
      delta_n = 0.1 * fmin(radi, radj);
    }
    return M_PI * delta_n * reff;
  }
}

/* ---------------------------------------------------------------------- */

double FixBatteryEIS::calculate_current_Li_SE(int i_Li, int j_SE, double phi_ed, double phi_el, double sigma_m)
{
  // Get equilibrium potential and exchange current density for Li particle
  double U_eq = 0.0;
  double i_0 = 0.01; // A/m2, placeholder value

  // Calculate overpotential for NMC cathode
  // η = φ_ed - φ_el - U_eq
  double eta = phi_ed - phi_el - U_eq;

  // Butler-Volmer equation with stress effect
  double RT = R * T;
  
  // Calculate exponential arguments and limit them
  double arg1 = (1.0 - alpha_a) * F * eta / RT - sigma_m * 9e-6 / RT;
  double arg2 = -alpha_c * F * eta / RT - sigma_m * 9e-6 / RT;
  
  // Prevent exponential overflow
  // if (arg1 > MAX_EXP_ARG) arg1 = MAX_EXP_ARG;
  // if (arg1 < -MAX_EXP_ARG) arg1 = -MAX_EXP_ARG;
  // if (arg2 > MAX_EXP_ARG) arg2 = MAX_EXP_ARG;
  // if (arg2 < -MAX_EXP_ARG) arg2 = -MAX_EXP_ARG;
  
  double exp_term1 = exp(arg1);
  double exp_term2 = exp(arg2);
  
  // Current density from Li to SE
  double i_pq = i_0 * (exp_term1 - exp_term2);
  
  // Limit current to prevent numerical issues
  // double i_max = 1000.0;  // A/m^2
  // if (i_pq > i_max) i_pq = i_max;
  // if (i_pq < -i_max) i_pq = -i_max;
  
  return i_pq;
}

/* ---------------------------------------------------------------------- */

double FixBatteryEIS::check_convergence()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  double local_error = 0.0;
  int count = 0;
  
  // Check convergence for both potentials
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      // Check electrolyte potential convergence for SE particles
      if (type[i] == SE_type) {
        double diff_el = fabs(phi_el[i] - phi_el_old[i]);
        if (diff_el > local_error) local_error = diff_el;
      }
      
      // Check electronic potential convergence for Li, SE, and CC particles
      // if (type[i] == Li_type || type[i] == SE_type || type[i] == BC_bottom_type) {
      //   double diff_ed = fabs(phi_ed[i] - phi_ed_old[i]);
      //   if (diff_ed > local_error) local_error = diff_ed;
      // }
      count++;
    }
  }
  
  double global_error;
  MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_MAX, world);
  
  return global_error;
}

/* ---------------------------------------------------------------------- */

double FixBatteryEIS::compute_scalar()
{
  // Return convergence error
  return convergence_error;
}

/* ---------------------------------------------------------------------- */

double FixBatteryEIS::compute_vector(int n)
{
  if (n == 0) return (double)current_iteration;
  else if (n == 1) return convergence_error;
  return 0.0;
}