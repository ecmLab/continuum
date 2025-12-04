/* ----------------------------------------------------------------------
    LIGGGHTS® DEM simulation engine
    
    Updates:
    - Calculate the current in (A) using the wanted current desnity and the cross-area of the xy plane.

    Notes: 
    - FOR CONTROLLED CURRENT in SYMMETRICAL CELL
    - Solves electrolyte potential with constant current BC at interfaces
    - Structure: CC_Atype / Li_Atype / SE / Li_Ctype / CC_Ctype
    - CORRECTED: Current conservation through cell
      - Total current is controlled at anode/SE interface
      - Same total current flows through SE/cathode interface
      - Current density at each interface depends on local contact area
------------------------------------------------------------------------- */

#include "fix_battery_eis.h"
#include "atom.h"
#include "update.h"
#include "compute.h"
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
#define MAX_EXP_ARG 20.0

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

  // Register property/atom for electronic potential (kept for compatibility)
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

  // Register property/atom for old electronic potential
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

  // Register property/atom for current density at interfaces
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
    fixarg[8]="0.0";
    fix_init_flag = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }
  
  // Register property/atom for per-particle contact area at anode interface
  fix_contact_area_anode = static_cast<FixPropertyAtom*>(modify->find_fix_property("contactAreaAnode","property/atom","scalar",0,0,style,false));
  if(!fix_contact_area_anode) {
    const char* fixarg[10];
    fixarg[0]="contactAreaAnode";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="contactAreaAnode";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_contact_area_anode = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }
  
  // Register property/atom for per-particle contact area at cathode interface
  fix_contact_area_cathode = static_cast<FixPropertyAtom*>(modify->find_fix_property("contactAreaCathode","property/atom","scalar",0,0,style,false));
  if(!fix_contact_area_cathode) {
    const char* fixarg[10];
    fixarg[0]="contactAreaCathode";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="contactAreaCathode";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_contact_area_cathode = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
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
  R(8.31),
  T(303.0),
  F(96485.0),
  sigma_el(0.05),
  sigma_Li(1.0e7),
  SE_type(1),
  Li_Ctype(2),
  Li_Atype(3),
  CC_Ctype(4),
  CC_Atype(5),
  phi_ref(0.0),
  cur_app(1.0),
  // New: Store calculated interface values
  total_current(0.0),
  global_area_anode(0.0),
  global_area_cathode(0.0),
  i_density_anode(0.0),
  i_density_cathode(0.0),
  phi_el(NULL),
  phi_el_old(NULL),
  phi_ed(NULL),
  phi_ed_old(NULL),
  current_Li_SE(NULL),
  hydrostatic_stress(NULL),
  contact_area_anode(NULL),
  contact_area_cathode(NULL),
  fix_phi_el(NULL),
  fix_phi_el_old(NULL),
  fix_phi_ed(NULL),
  fix_phi_ed_old(NULL),
  fix_current_Li_SE(NULL),
  fix_hydrostatic_stress(NULL),
  fix_init_flag(NULL),
  fix_contact_area_anode(NULL),
  fix_contact_area_cathode(NULL),
  list(NULL),
  first_run(true)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix battery/eis command");

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
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix battery/eis command - BC_types needs 4 values");
      Li_Ctype = force->inumeric(FLERR,arg[iarg+1]);
      Li_Atype = force->inumeric(FLERR,arg[iarg+2]);
      CC_Ctype = force->inumeric(FLERR,arg[iarg+3]);
      CC_Atype = force->inumeric(FLERR,arg[iarg+4]);
      iarg += 5;
    } else if (strcmp(arg[iarg],"BC_potentials") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix battery/eis command");
      phi_ref = force->numeric(FLERR,arg[iarg+1]);
      cur_app = force->numeric(FLERR,arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"conductivity") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix battery/eis command");
      sigma_el = force->numeric(FLERR,arg[iarg+1]);
      sigma_Li = force->numeric(FLERR,arg[iarg+2]);
      iarg += 3;
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
  size_vector = 6;  // Extended to include interface info
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
  if (force->newton_pair == 1)
    error->all(FLERR,"Fix battery/eis requires newton pair off");

  fix_phi_el = static_cast<FixPropertyAtom*>(modify->find_fix_property("electrolytePotential","property/atom","scalar",0,0,style));
  fix_phi_el_old = static_cast<FixPropertyAtom*>(modify->find_fix_property("electrolytePotentialOld","property/atom","scalar",0,0,style));
  fix_phi_ed = static_cast<FixPropertyAtom*>(modify->find_fix_property("electronicPotential","property/atom","scalar",0,0,style));
  fix_phi_ed_old = static_cast<FixPropertyAtom*>(modify->find_fix_property("electronicPotentialOld","property/atom","scalar",0,0,style));
  fix_current_Li_SE = static_cast<FixPropertyAtom*>(modify->find_fix_property("currentLiSE","property/atom","scalar",0,0,style));
  fix_hydrostatic_stress = static_cast<FixPropertyAtom*>(modify->find_fix_property("hydrostaticStress","property/atom","scalar",0,0,style));
  fix_init_flag = static_cast<FixPropertyAtom*>(modify->find_fix_property("batteryInitFlag","property/atom","scalar",0,0,style));
  fix_contact_area_anode = static_cast<FixPropertyAtom*>(modify->find_fix_property("contactAreaAnode","property/atom","scalar",0,0,style));
  fix_contact_area_cathode = static_cast<FixPropertyAtom*>(modify->find_fix_property("contactAreaCathode","property/atom","scalar",0,0,style));

  if(!fix_phi_el || !fix_phi_el_old || !fix_phi_ed || !fix_phi_ed_old || 
     !fix_current_Li_SE || !fix_hydrostatic_stress || !fix_init_flag ||
     !fix_contact_area_anode || !fix_contact_area_cathode)
    error->all(FLERR,"Could not find required property/atom fixes");

  if (SE_type < 1 || SE_type > atom->ntypes)
    error->all(FLERR,"Invalid SE particle type for battery/eis");
  if (Li_Ctype < 1 || Li_Ctype > atom->ntypes || Li_Atype < 1 || Li_Atype > atom->ntypes)
    error->all(FLERR,"Invalid Li particle types for battery/eis");
  if (CC_Ctype < 1 || CC_Ctype > atom->ntypes || CC_Atype < 1 || CC_Atype > atom->ntypes)
    error->all(FLERR,"Invalid CC particle types for battery/eis");

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

  bool already_initialized = false;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && init_flag[i] > 0.5) {
      already_initialized = true;
      break;
    }
  }
  
  if (!already_initialized) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (type[i] == SE_type) {
          phi_el[i] = phi_ref;
          phi_el_old[i] = phi_ref;
        } else if (type[i] == Li_Ctype || type[i] == Li_Atype) {
          phi_el[i] = phi_ref;
          phi_el_old[i] = phi_ref;
        } else if (type[i] == CC_Ctype) {
          phi_el[i] = phi_ref;
          phi_el_old[i] = phi_ref;
        } else if (type[i] == CC_Atype) {
          phi_el[i] = phi_ref;
          phi_el_old[i] = phi_ref;
        }
        phi_ed[i] = 0.0;
        phi_ed_old[i] = 0.0;
        init_flag[i] = 1.0;
      }
    }
    first_run = false;
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        phi_el_old[i] = phi_el[i];
        phi_ed_old[i] = phi_ed[i];
      }
    }
  }
  
  apply_reference_potential();
}

/* ---------------------------------------------------------------------- */

void FixBatteryEIS::pre_force(int vflag)
{
  updatePtrs();
  update->vflag_atom = update->ntimestep;
  calculate_hydrostatic_stress();
}

/* ---------------------------------------------------------------------- */

void FixBatteryEIS::post_force(int vflag)
{
  // First, calculate interface contact areas and current distribution
  calculate_interface_currents();
  
  // Then run SOR iterations
  for (int iter = 0; iter < max_iterations; iter++) {
    solve_eis_iteration();
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
  contact_area_anode = fix_contact_area_anode->vector_atom;
  contact_area_cathode = fix_contact_area_cathode->vector_atom;
}

/* ---------------------------------------------------------------------- */
/* NEW FUNCTION: Calculate interface areas and current distribution
   
   Key physics:
   1. Applied current density (cur_app) is specified at anode/SE interface
   2. Total current I_total = cur_app * A_anode_SE
   3. This same total current must flow through cathode/SE interface
   4. Current density at cathode = I_total / A_cathode_SE
   
   This ensures current conservation: what goes in must come out
------------------------------------------------------------------------- */

void FixBatteryEIS::calculate_interface_currents()
{
  int i,j,ii,jj,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,r;
  
  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // Reset per-particle contact areas
  for (i = 0; i < nlocal; i++) {
    contact_area_anode[i] = 0.0;
    contact_area_cathode[i] = 0.0;
  }

  // Local accumulators for total interface areas
  double local_area_anode_SE = 0.0;   // Li_Atype to SE_type interface
  double local_area_cathode_SE = 0.0; // Li_Ctype to SE_type interface
  
  // First pass: Calculate per-particle and total contact areas at each interface
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      r = sqrt(rsq);
      
      if (r < (radius[i] + radius[j])) {
        double contact_area = calculate_contact_area(i, j);
        
        if (contact_area > 0.0) {
          // Anode interface: Li_Atype <-> SE_type
          if (type[i] == Li_Atype && type[j] == SE_type) {
            contact_area_anode[i] += contact_area;
            local_area_anode_SE += contact_area;
          }
          else if (type[i] == SE_type && type[j] == Li_Atype) {
            contact_area_anode[i] += contact_area;
            // Don't double count in total - only count from Li side
          }
          
          // Cathode interface: Li_Ctype <-> SE_type
          if (type[i] == Li_Ctype && type[j] == SE_type) {
            contact_area_cathode[i] += contact_area;
            local_area_cathode_SE += contact_area;
          }
          else if (type[i] == SE_type && type[j] == Li_Ctype) {
            contact_area_cathode[i] += contact_area;
            // Don't double count in total - only count from Li side
          }
        }
      }
    }
  }
  
  // Global reduction for total contact areas
  MPI_Allreduce(&local_area_anode_SE, &global_area_anode, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&local_area_cathode_SE, &global_area_cathode, 1, MPI_DOUBLE, MPI_SUM, world);
  
  // Calculate cycling direction based on timestep
  double ndt = update->ntimestep;
  int period = static_cast<int>(ndt) / 36000; // Assuming 36000 timesteps for 60 mins (60 min charge/discharge)
  double sign = (period % 2 == 0) ? 1.0 : -1.0;

  // CCD Test Protocol: Increase current magnitude every two periods
  int num_positive_steps = period / 2;
  double increase_step = cur_app / 4.0;
  double current_magnitude = cur_app + (num_positive_steps * increase_step);
  double i_density = sign * current_magnitude; // Applied current density in A/m2 // sign * current_magnitude
  
  // // Calculate applied current density
  // int period = static_cast<int>(ndt) / 3600; // Assuming 3600 is 1 hour
  // double sign = (period % 2 == 0) ? 1.0 : -1.0;
  // int num_positive_steps = period / 2;
  // double increase_step = cur_app / 4.0;
  // double current_magnitude = cur_app + (num_positive_steps * increase_step);
  // double i_app = current_magnitude; // Applied current density in A/m2 // sign * current_magnitude


  // Applied current density at anode (controlled)
  // double i_density = cur_app;  // A/m² sign * cur_app

  total_current = i_density * 4e-10; // A Ideal current assuming area is 20 um x 20 um
  
  // Calculate total current from anode interface
  // i_density_anode = I_total / A_anode
  if (global_area_anode > SMALL) {
    i_density_anode = total_current / global_area_anode;  // Amperes
  } else {
    i_density_anode = 0.0;
  }
  
  // Calculate current density at cathode interface
  // i_density_cathode = I_total / A_cathode
  // This ensures current conservation through the cell
  if (global_area_cathode > SMALL) {
    i_density_cathode = total_current / global_area_cathode;  // A/m²
  } else {
    i_density_cathode = 0.0;
  }
  
  // Debug output (optional - can be removed in production)
  if (update->ntimestep % 720 == 0) {
    if (comm->me == 0) {
      fprintf(screen, "Battery EIS - Step %ld:\n", update->ntimestep);
      fprintf(screen, "  Anode interface:   A = %.6e m², i = %.4f A/m²\n", 
              global_area_anode, i_density_anode);
      fprintf(screen, "  Cathode interface: A = %.6e m², i = %.4f A/m²\n", 
              global_area_cathode, i_density_cathode);
      fprintf(screen, "  Total current:     I = %.6e A\n", total_current);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixBatteryEIS::solve_eis_iteration()
{
  int i,j,ii,jj,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,r,r_SI;
  double phi_sum,coeff_sum,cur_source;
  
  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // Store old values for convergence check
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      phi_el_old[i] = phi_el[i];
    }
  }

  // Reset current tracking
  for (i = 0; i < nlocal; i++) {
    current_Li_SE[i] = 0.0;
  }

  // Main SOR iteration loop
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    
    // Skip CC_Ctype particles - these are the reference electrode (ground)
    // if (type[i] == CC_Ctype) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    
    phi_sum = 0.0;
    coeff_sum = 0.0;
    cur_source = 0.0;
    
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      r = sqrt(rsq);
      r_SI = r * 1.0e-6;
      
      if (r < (radius[i] + radius[j])) {
        double contact_area = calculate_contact_area(i, j);
        
        if (contact_area > 0.0) {
          double sigma_eff = 0.0;
          
          // Determine effective conductivity based on interface type
          if (type[i] == SE_type && type[j] == SE_type) {
            sigma_eff = sigma_el;
          }
          else if ((type[i] == Li_Atype && type[j] == Li_Atype) ||
                   (type[i] == Li_Ctype && type[j] == Li_Ctype)) {
            sigma_eff = sigma_Li;
          }
          else if ((type[i] == SE_type && (type[j] == Li_Atype || type[j] == Li_Ctype)) ||
                   ((type[i] == Li_Atype || type[i] == Li_Ctype) && type[j] == SE_type)) {
            sigma_eff = 2.0 * sigma_el * sigma_Li / (sigma_el + sigma_Li);
          }
          else if ((type[i] == CC_Atype && type[j] == Li_Atype) ||
                   (type[i] == Li_Atype && type[j] == CC_Atype) ||
                   (type[i] == CC_Ctype && type[j] == Li_Ctype) ||
                   (type[i] == Li_Ctype && type[j] == CC_Ctype)) {
            sigma_eff = sigma_Li;
          }
          else if ((type[i] == CC_Atype && type[j] == CC_Atype) ||
                   (type[i] == CC_Ctype && type[j] == CC_Ctype)) {
            sigma_eff = sigma_Li * 100.0;
          }
          
          if (sigma_eff > SMALL) {
            double conductance = sigma_eff * contact_area / r_SI;
            phi_sum += conductance * phi_el[j];
            coeff_sum += conductance;
          }
          
          /* CORRECTED CURRENT SOURCE TERMS
           *
           * At anode interface (Li_Atype <-> SE_type):
           *   - Apply controlled current density i_density_anode
           *   - Current enters the electrolyte from lithium (positive source for Li)
           *
           * At cathode interface (Li_Ctype <-> SE_type):
           *   - Apply calculated current density i_density_cathode
           *   - This ensures total current conservation
           *   - Current exits the electrolyte to lithium (negative source for Li)
           *
           * The key insight: we don't apply the same current DENSITY at both
           * interfaces. We apply the same total CURRENT, distributed according
           * to the local contact area.
           */
          
          // Anode interface: Li_Atype particle contacting SE_type
          if (type[i] == Li_Atype && type[j] == SE_type) {
            // Current injection at anode (positive = stripping Li)
            double local_current = i_density_anode * contact_area;
            cur_source += local_current;
            current_Li_SE[i] += local_current;
          }
          // SE particle at anode interface (receives current from Li_Atype)
          // else if (type[i] == SE_type && type[j] == Li_Atype) {
          //   double local_current = -1 * i_density_anode * contact_area;
          //   cur_source += local_current;  // Current flows into SE
          //   current_Li_SE[i] += local_current;
          // }
          
          // Cathode interface: Li_Ctype particle contacting SE_type
          else if (type[i] == Li_Ctype && type[j] == SE_type) {
            // Current extraction at cathode (negative = plating Li)
            double local_current = -1 * i_density_cathode * contact_area;
            cur_source += local_current;
            current_Li_SE[i] += local_current;
          }
          // SE particle at cathode interface (releases current to Li_Ctype)
          // else if (type[i] == SE_type && type[j] == Li_Ctype) {
          //   double local_current = i_density_cathode * contact_area;
          //   cur_source += local_current;  // Current flows out of SE
          //   current_Li_SE[i] += local_current;
          // }
        }
      }
    }
    
    if (coeff_sum > SMALL) {
      double phi_new = (phi_sum + cur_source) / coeff_sum;
      
      if (!std::isnan(phi_new) && !std::isinf(phi_new)) {
        phi_el[i] = phi_el_old[i] + omega * (phi_new - phi_el_old[i]);
      } else {
        phi_el[i] = phi_el_old[i];
      }
    }
  }
  
  // Forward communication
  fix_phi_el->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

void FixBatteryEIS::apply_reference_potential()
{
  // Reference potential at bottom current collector (ground)
  // All other potentials evolve based on current flow
  
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (type[i] == CC_Ctype) {
        phi_el[i] = phi_ref;
        phi_el_old[i] = phi_ref;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixBatteryEIS::calculate_hydrostatic_stress()
{
  char *stress_id = (char *)"st"; 
  int icompute = modify->find_compute(stress_id);
  
  if (icompute < 0) {
    error->all(FLERR, "FixBatteryEIS: Could not find compute with ID 'st'.");
  }
  
  Compute *stress_compute = modify->compute[icompute];

  if (!(stress_compute->invoked_flag & INVOKED_PERATOM)) {
      stress_compute->compute_peratom();
      stress_compute->invoked_flag |= INVOKED_PERATOM;
  }

  double **stress = stress_compute->array_atom;

  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  double stress_conversion = 1.0e-15;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && (type[i] == Li_Atype || type[i] == Li_Ctype)) {
      double radius_SI = radius[i] * 1.0e-6;
      double volume_SI = (4.0 / 3.0) * M_PI * radius_SI * radius_SI * radius_SI;
      double trace_virial_micro = stress[i][0] + stress[i][1] + stress[i][2];
      double trace_virial_SI = trace_virial_micro * stress_conversion;

      if (volume_SI > 0.0) {
        hydrostatic_stress[i] = trace_virial_SI / (3.0 * volume_SI); // trace_virial_SI / (3.0 * volume_SI)
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
  double r = sqrt(delx*delx + dely*dely + delz*delz);
  
  double radi = radius[i] * length_conversion;
  double radj = radius[j] * length_conversion;
  double radsum = radi + radj;
  
  if (r >= radsum) return 0.0;
  
  if (r < fmax(radi, radj)) {
    double rmin = fmin(radi, radj);
    return M_PI * rmin * rmin;
  } else {
    double delta_n = radsum - r;
    double reff = radi * radj / radsum;
    if (delta_n > 0.1 * fmin(radi, radj)) {
      delta_n = 0.1 * fmin(radi, radj);
    }
    return M_PI * delta_n * reff;
  }
}

/* ---------------------------------------------------------------------- */

double FixBatteryEIS::calculate_current_Li_SE(int i_Li, int j_SE, double phi_ed, double phi_el, double sigma_m)
{
  // Not used in constant current mode, kept for compatibility
  double U_eq = 0.0;
  double i_0 = 0.01;
  double eta = phi_ed - phi_el - U_eq;
  double RT = R * T;
  
  double arg1 = (1.0 - alpha_a) * F * eta / RT - sigma_m * 9e-6 / RT;
  double arg2 = -alpha_c * F * eta / RT - sigma_m * 9e-6 / RT;
  
  double exp_term1 = exp(arg1);
  double exp_term2 = exp(arg2);
  
  return i_0 * (exp_term1 - exp_term2);
}

/* ---------------------------------------------------------------------- */

double FixBatteryEIS::check_convergence()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  double local_error = 0.0;
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (type[i] == SE_type || type[i] == Li_Atype || type[i] == Li_Ctype) {
        double diff_el = fabs(phi_el[i] - phi_el_old[i]);
        if (diff_el > local_error) local_error = diff_el;
      }
    }
  }
  
  double global_error;
  MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_MAX, world);
  
  return global_error;
}

/* ---------------------------------------------------------------------- */

double FixBatteryEIS::compute_scalar()
{
  return convergence_error;
}

/* ---------------------------------------------------------------------- */

double FixBatteryEIS::compute_vector(int n)
{
  // Extended output for monitoring
  if (n == 0) return (double)current_iteration;
  else if (n == 1) return convergence_error;
  else if (n == 2) return global_area_anode;      // Anode interface area
  else if (n == 3) return global_area_cathode;    // Cathode interface area
  else if (n == 4) return i_density_anode;        // Current density at anode
  else if (n == 5) return i_density_cathode;      // Current density at cathode
  return 0.0;
}