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
    
    Copyright 2024-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include "fix_battery_sor.h"
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

#define SMALL 1e-10
#define MAX_EXP_ARG 20.0  // Prevent exp() overflow

/* ---------------------------------------------------------------------- */

void FixBatterySOR::post_create()
{
  // Register property/atom for electric potential
  fix_phi_el = static_cast<FixPropertyAtom*>(modify->find_fix_property("electricPotential","property/atom","scalar",0,0,style,false));
  if(!fix_phi_el) {
    const char* fixarg[10];
    fixarg[0]="electricPotential";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="electricPotential";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_phi_el = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  // Register property/atom for old electric potential (convergence check)
  fix_phi_el_old = static_cast<FixPropertyAtom*>(modify->find_fix_property("electricPotentialOld","property/atom","scalar",0,0,style,false));
  if(!fix_phi_el_old) {
    const char* fixarg[10];
    fixarg[0]="electricPotentialOld";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="electricPotentialOld";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_phi_el_old = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  // Register property/atom for current density AM-SE
  fix_current_AM_SE = static_cast<FixPropertyAtom*>(modify->find_fix_property("currentAMSE","property/atom","scalar",0,0,style,false));
  if(!fix_current_AM_SE) {
    const char* fixarg[10];
    fixarg[0]="currentAMSE";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="currentAMSE";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_current_AM_SE = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
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
}

/* ---------------------------------------------------------------------- */

FixBatterySOR::FixBatterySOR(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  omega(1.9),
  tolerance(1e-6),
  max_iterations(10000),
  current_iteration(0),
  convergence_error(0.0),
  R(8.31), // J*mol^(-1)*K^(-1)
  T(303.0), // K
  F(96485.0), // C/mol
  sigma_SE(0.05), // S/m
  alpha_a(0.5), //unitless
  alpha_c(0.5),
  BC_bottom_type(2),      // Default: CC1
  BC_top_type(6),         // Default: CC2
  phi_el_BC_bottom(0.0),  // Default: 0V at bottom
  phi_el_BC_top(0.04),    // Default: 0.04V at top
  phi_el(NULL),
  phi_el_old(NULL),
  equilibrium_potential(NULL),
  exchange_current_density(NULL),
  current_AM_SE(NULL),
  hydrostatic_stress(NULL),
  fix_phi_el(NULL),
  fix_phi_el_old(NULL),
  fix_equilibrium_potential(NULL),
  fix_exchange_current_density(NULL),
  fix_current_AM_SE(NULL),
  fix_hydrostatic_stress(NULL),
  groupbit_SE(0),
  groupbit_AM(0),
  SE_type(3),
  AM_type(1),
  list(NULL)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix battery/sor command");

  // Parse arguments
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"omega") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix battery/sor command");
      omega = force->numeric(FLERR,arg[iarg+1]);
      if (omega <= 0.0 || omega >= 2.0)
        error->all(FLERR,"SOR omega must be between 0 and 2");
      iarg += 2;
    } else if (strcmp(arg[iarg],"tolerance") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix battery/sor command");
      tolerance = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"max_iter") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix battery/sor command");
      max_iterations = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"temperature") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix battery/sor command");
      T = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"BC_types") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix battery/sor command");
      BC_bottom_type = force->inumeric(FLERR,arg[iarg+1]);
      BC_top_type = force->inumeric(FLERR,arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"BC_potentials") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix battery/sor command");
      phi_el_BC_bottom = force->numeric(FLERR,arg[iarg+1]);
      phi_el_BC_top = force->numeric(FLERR,arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"SE_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix battery/sor command");
      SE_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"AM_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix battery/sor command");
      AM_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix battery/sor command");
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

FixBatterySOR::~FixBatterySOR()
{
}

/* ---------------------------------------------------------------------- */

int FixBatterySOR::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBatterySOR::init()
{
  // Check newton pair
  if (force->newton_pair == 1)
    error->all(FLERR,"Fix battery/sor requires newton pair off");

  // Find required fixes
  fix_equilibrium_potential = static_cast<FixPropertyAtom*>(modify->find_fix_property("equilibriumPotential","property/atom","scalar",0,0,style));
  if(!fix_equilibrium_potential)
    error->all(FLERR,"Fix battery/sor requires equilibriumPotential property");
    
  fix_exchange_current_density = static_cast<FixPropertyAtom*>(modify->find_fix_property("exchangeCurrentDensity","property/atom","scalar",0,0,style));
  if(!fix_exchange_current_density)
    error->all(FLERR,"Fix battery/sor requires exchangeCurrentDensity property");

  // Get all property fixes
  fix_phi_el = static_cast<FixPropertyAtom*>(modify->find_fix_property("electricPotential","property/atom","scalar",0,0,style));
  fix_phi_el_old = static_cast<FixPropertyAtom*>(modify->find_fix_property("electricPotentialOld","property/atom","scalar",0,0,style));
  fix_current_AM_SE = static_cast<FixPropertyAtom*>(modify->find_fix_property("currentAMSE","property/atom","scalar",0,0,style));
  fix_hydrostatic_stress = static_cast<FixPropertyAtom*>(modify->find_fix_property("hydrostaticStress","property/atom","scalar",0,0,style));

  if(!fix_phi_el || !fix_phi_el_old || !fix_current_AM_SE || !fix_hydrostatic_stress)
    error->all(FLERR,"Could not find required property/atom fixes");

  // Validate particle types
  if (SE_type < 1 || SE_type > atom->ntypes || AM_type < 1 || AM_type > atom->ntypes)
    error->all(FLERR,"Invalid particle types for battery/sor");
    
  if (BC_bottom_type < 1 || BC_bottom_type > atom->ntypes || BC_top_type < 1 || BC_top_type > atom->ntypes)
    error->all(FLERR,"Invalid boundary condition types for battery/sor");

  // Request neighbor list
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixBatterySOR::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixBatterySOR::setup(int vflag)
{
  updatePtrs();
  
  int nlocal = atom->nlocal;
  int *type = atom->type;
  int *mask = atom->mask;
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      // Initialize based on particle type
      if (type[i] == SE_type) {
        phi_el[i] = 0.0;
        phi_el_old[i] = 0.0;
      } else if (type[i] == AM_type) {
        // AM particles don't have a direct electric potential
        phi_el[i] = 0.0;
        phi_el_old[i] = 0.0;
      } else if (type[i] == BC_bottom_type) {
        // Bottom boundary condition (e.g., CC1)
        phi_el[i] = phi_el_BC_bottom;
        phi_el_old[i] = phi_el_BC_bottom;
      } else if (type[i] == BC_top_type) {
        // Top boundary condition (e.g., CC2)
        phi_el[i] = phi_el_BC_top;
        phi_el_old[i] = phi_el_BC_top;
      }
    }
  }
  
  // Apply boundary conditions
  apply_boundary_conditions();
}

/* ---------------------------------------------------------------------- */

void FixBatterySOR::pre_force(int vflag)
{
  updatePtrs();
  
  // Calculate hydrostatic stress for AM particles
  calculate_hydrostatic_stress();
}

/* ---------------------------------------------------------------------- */

void FixBatterySOR::post_force(int vflag)
{
  // Solve using SOR iterations
  current_iteration = 0;
  convergence_error = 1.0;
  
  while (current_iteration < max_iterations && convergence_error > tolerance) {
    solve_sor_iteration();
    convergence_error = check_convergence();
    current_iteration++;
  }
  
  if (current_iteration >= max_iterations && comm->me == 0) {
    error->warning(FLERR,"Battery SOR did not converge");
  }
}

/* ---------------------------------------------------------------------- */

void FixBatterySOR::updatePtrs()
{
  phi_el = fix_phi_el->vector_atom;
  phi_el_old = fix_phi_el_old->vector_atom;
  equilibrium_potential = fix_equilibrium_potential->vector_atom;
  exchange_current_density = fix_exchange_current_density->vector_atom;
  current_AM_SE = fix_current_AM_SE->vector_atom;
  hydrostatic_stress = fix_hydrostatic_stress->vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixBatterySOR::solve_sor_iteration()
{
  int i,j,ii,jj,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,r,r_SI;
  double phi_sum,coeff_sum,cur_sum;
  
  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // Store old values for SE particles only
  for (i = 0; i < nlocal; i++) {
    if (type[i] == SE_type) {
      phi_el_old[i] = phi_el[i];
    }
  }

  // Reset current_AM_SE for all particles
  for (int i = 0; i < nlocal; i++) {
    current_AM_SE[i] = 0.0;
  }

  // SOR iteration - ONLY for SE particles
  for (ii = 0; ii < inum; ii++) { // 1
    i = ilist[ii];
    
    // Only update SE particles, skip AM and boundary condition particles
    if (type[i] != SE_type) continue;
    
    // For SE particles: solve diffusion equation
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    
    phi_sum = 0.0;
    coeff_sum = 0.0;
    cur_sum = 0.0;
    
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    for (jj = 0; jj < jnum; jj++) { // 2
      j = jlist[jj];
      j &= NEIGHMASK;
      
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      r = sqrt(rsq);  // micrometers;
      r_SI = r * 1.0e-6;  // Convert to m;      

      if (r < (radius[i] + radius[j])) { // 3
        double contact_area = calculate_contact_area(i, j); //Should be m2
        
        if (contact_area > 0.0) { // 4
          if (type[j] == SE_type) {
            // SE-SE conductivity
            double conductance = sigma_SE * contact_area / r_SI;
            phi_sum += conductance * phi_el[j];
            coeff_sum += conductance;
            
          } else if (type[j] == AM_type) {
            // SE-AM interface: calculate current from AM to SE
            double phi_AM_effective = equilibrium_potential[j];
            
            // Calculate current using Butler-Volmer equation
            double i_pq = calculate_current_AM_SE(j, i, phi_AM_effective, phi_el[i], hydrostatic_stress[j]); // Current density A/m2
            current_AM_SE[j] += i_pq * contact_area; // Local current in a AM particle due to electro-chemical reaction
            
            // Add current contribution to SE particle
            cur_sum += i_pq * contact_area; // Local current in a SE particle due to electro-chemical reaction
            
          } else if (type[j] == BC_bottom_type || type[j] == BC_top_type) {
            // SE in contact with boundary condition particle (CC1 or CC2)
            // These act as Dirichlet boundary conditions
            double conductance = sigma_SE * contact_area / r_SI;
            phi_sum += conductance * phi_el[j];  // phi_el[j] is fixed at BC value
            coeff_sum += conductance;
          }
        } // 4
      } // 3
    } // 2
    
    // Update SE potential using SOR
    if (coeff_sum > SMALL) {
      double phi_new = (phi_sum + cur_sum) / coeff_sum; // Updated 7/30/2025
      
      // Check for NaN
      if (std::isnan(phi_new)) { // Old: phi_new != phi_new
        error->warning(FLERR,"NaN detected in SE potential calculation");
        phi_el[i] = phi_el_old[i];
      } else {
        phi_el[i] = phi_el_old[i] + omega * (phi_new - phi_el_old[i]);
      }
    }
  } // 1
  
  // Apply boundary conditions after each iteration
  apply_boundary_conditions();
  
  // Forward communication for SE potentials
  fix_phi_el->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

void FixBatterySOR::apply_boundary_conditions()
{
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      // Bottom boundary condition particles (e.g., CC1)
      if (type[i] == BC_bottom_type) {
        phi_el[i] = phi_el_BC_bottom;
        phi_el_old[i] = phi_el_BC_bottom;
      }
      // Top boundary condition particles (e.g., CC2)
      else if (type[i] == BC_top_type) {
        phi_el[i] = phi_el_BC_top;
        phi_el_old[i] = phi_el_BC_top;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixBatterySOR::calculate_hydrostatic_stress()
{
  double **f = atom->f;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && type[i] == AM_type) {
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

double FixBatterySOR::calculate_contact_area(int i, int j)
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

double FixBatterySOR::calculate_current_AM_SE(int i_AM, int j_SE, double phi_el_AM, double phi_el_SE, double sigma_m)
{
  // Get equilibrium potential and exchange current density for AM particle
  double U_eq = equilibrium_potential[i_AM];
  double i_0 = exchange_current_density[i_AM];
  
  // Calculate overpotential
  // Note: phi_el_AM should be the equilibrium potential for AM particles
  double eta = 0 - phi_el_SE - U_eq;
  
  // Limit overpotential to prevent numerical overflow
  // if (eta > 0.9) eta = 0.9;
  // if (eta < -0.9) eta = -0.9;
  
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
  
  // Current density from AM to SE
  double i_pq = i_0 * (exp_term1 - exp_term2);
  
  // Limit current to prevent numerical issues
  // double i_max = 1000.0;  // A/m^2
  // if (i_pq > i_max) i_pq = i_max;
  // if (i_pq < -i_max) i_pq = -i_max;
  
  return i_pq;
}

/* ---------------------------------------------------------------------- */

double FixBatterySOR::check_convergence()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  double local_error = 0.0;
  int count = 0;
  
  // Only check convergence for SE particles
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && type[i] == SE_type) {
      double diff = fabs(phi_el[i] - phi_el_old[i]);
      if (diff > local_error) local_error = diff;
      count++;
    }
  }
  
  double global_error;
  MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_MAX, world);
  
  return global_error;
}

/* ---------------------------------------------------------------------- */

double FixBatterySOR::compute_scalar()
{
  // Return convergence error
  return convergence_error;
}

/* ---------------------------------------------------------------------- */

double FixBatterySOR::compute_vector(int n)
{
  if (n == 0) return (double)current_iteration;
  else if (n == 1) return convergence_error;
  return 0.0;
}
