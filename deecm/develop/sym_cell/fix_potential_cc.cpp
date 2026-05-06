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
    - FOR CONTROLLED CURRENT
    - Solves both electronic and electrolyte potentials
------------------------------------------------------------------------- */

#include "fix_potential_cc.h"
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

void FixPotentialCC::post_create()
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
  fix_init_flag = static_cast<FixPropertyAtom*>(modify->find_fix_property("potentialInitFlag","property/atom","scalar",0,0,style,false));
  if(!fix_init_flag) {
    const char* fixarg[10];
    fixarg[0]="potentialInitFlag";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="potentialInitFlag";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.0";
    fix_init_flag = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }
}

/* ---------------------------------------------------------------------- */

FixPotentialCC::FixPotentialCC(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  omega(1.9),
  tolerance(1e-6),
  convergence_error(0.0),
  R(8.31),                   // J*mol^(-1)*K^(-1)  (constant, not user-configurable)
  T(303.0),                  // K
  F(96485.0),                // C/mol              (constant, not user-configurable)
  sigma_el(0.05),            // S/m  - Electrolyte ionic conductivity
  sigma_ed_SE(200.0),        // S/m  - Electronic conductivity for CBD/SE
  sigma_ed_CC(300.0),        // S/m  - Electronic conductivity for CC
  alpha_a(0.5),
  alpha_c(0.5),
  i_0(0.01),                 // A/m² - Exchange current density
  U_eq(0.0),                 // V    - Equilibrium potential
  Li_Ctype(2),
  Li_Atype(3),
  phi_el_BC_Cat(0.0),        // V
  phi_el_BC_An(0.01),        // V
  phi_ed_BC_anode(0.0),      // V
  cur_app(0.0),              // A/m²
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
    error->all(FLERR,"Illegal fix potential/cc command");

  // Parse arguments
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"omega") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix potential/cc command: omega requires 1 value");
      omega = force->numeric(FLERR,arg[iarg+1]);
      if (omega <= 0.0 || omega >= 2.0)
        error->all(FLERR,"CC omega must be between 0 and 2");
      iarg += 2;
    } else if (strcmp(arg[iarg],"tolerance") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix potential/cc command: tolerance requires 1 value");
      tolerance = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"temperature") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix potential/cc command: temperature requires 1 value");
      T = force->numeric(FLERR,arg[iarg+1]);
      if (T <= 0.0)
        error->all(FLERR,"Temperature must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg],"BC_types") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix potential/cc command: BC_types requires 2 values");
      Li_Ctype = force->inumeric(FLERR,arg[iarg+1]);
      Li_Atype = force->inumeric(FLERR,arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"BC_potentials") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix potential/cc command: BC_potentials requires 2 values");
      phi_el_BC_Cat = force->numeric(FLERR,arg[iarg+1]);
      phi_el_BC_An = force->numeric(FLERR,arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"conductivity") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix potential/cc command: conductivity requires 1 value");
      sigma_el = force->numeric(FLERR,arg[iarg+1]);
      if (sigma_el <= 0.0)
        error->all(FLERR,"Electrolyte conductivity must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg],"sigma_ed_SE") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix potential/cc command: sigma_ed_SE requires 1 value");
      sigma_ed_SE = force->numeric(FLERR,arg[iarg+1]);
      if (sigma_ed_SE <= 0.0)
        error->all(FLERR,"sigma_ed_SE must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg],"sigma_ed_CC") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix potential/cc command: sigma_ed_CC requires 1 value");
      sigma_ed_CC = force->numeric(FLERR,arg[iarg+1]);
      if (sigma_ed_CC <= 0.0)
        error->all(FLERR,"sigma_ed_CC must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg],"SE_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix potential/cc command: SE_type requires 1 value");
      SE_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"alpha") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix potential/cc command: alpha requires 2 values (alpha_a alpha_c)");
      alpha_a = force->numeric(FLERR,arg[iarg+1]);
      alpha_c = force->numeric(FLERR,arg[iarg+2]);
      if (alpha_a <= 0.0 || alpha_c <= 0.0)
        error->all(FLERR,"alpha_a and alpha_c must be positive");
      iarg += 3;
    } else if (strcmp(arg[iarg],"i_0") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix potential/cc command: i_0 requires 1 value");
      i_0 = force->numeric(FLERR,arg[iarg+1]);
      if (i_0 <= 0.0)
        error->all(FLERR,"i_0 must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg],"U_eq") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix potential/cc command: U_eq requires 1 value");
      U_eq = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"phi_ed_anode") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix potential/cc command: phi_ed_anode requires 1 value");
      phi_ed_BC_anode = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"cur_app") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix potential/cc command: cur_app requires 1 value");
      cur_app = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix potential/cc command: unknown keyword");
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

FixPotentialCC::~FixPotentialCC()
{
}

/* ---------------------------------------------------------------------- */

int FixPotentialCC::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPotentialCC::init()
{
  // Check newton pair
  if (force->newton_pair == 1)
    error->all(FLERR,"Fix potential/cc requires newton pair off");

  // Get all property fixes
  fix_phi_el = static_cast<FixPropertyAtom*>(modify->find_fix_property("electrolytePotential","property/atom","scalar",0,0,style));
  fix_phi_el_old = static_cast<FixPropertyAtom*>(modify->find_fix_property("electrolytePotentialOld","property/atom","scalar",0,0,style));
  fix_phi_ed = static_cast<FixPropertyAtom*>(modify->find_fix_property("electronicPotential","property/atom","scalar",0,0,style));
  fix_phi_ed_old = static_cast<FixPropertyAtom*>(modify->find_fix_property("electronicPotentialOld","property/atom","scalar",0,0,style));
  fix_current_Li_SE = static_cast<FixPropertyAtom*>(modify->find_fix_property("currentLiSE","property/atom","scalar",0,0,style));
  fix_hydrostatic_stress = static_cast<FixPropertyAtom*>(modify->find_fix_property("hydrostaticStress","property/atom","scalar",0,0,style));
  fix_init_flag = static_cast<FixPropertyAtom*>(modify->find_fix_property("potentialInitFlag","property/atom","scalar",0,0,style));

  if(!fix_phi_el || !fix_phi_el_old || !fix_phi_ed || !fix_phi_ed_old || !fix_current_Li_SE || !fix_hydrostatic_stress || !fix_init_flag)
    error->all(FLERR,"Could not find required property/atom fixes");

  // Validate particle types
  if (SE_type < 1 || SE_type > atom->ntypes)
    error->all(FLERR,"Invalid particle types for potential/cc");
    
  if (Li_Ctype < 1 || Li_Ctype > atom->ntypes || Li_Atype < 1 || Li_Atype > atom->ntypes)
    error->all(FLERR,"Invalid boundary condition types for potential/cc");

  // Request neighbor list
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixPotentialCC::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixPotentialCC::setup(int vflag)
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
        if (type[i] == SE_type) {
          phi_el[i] = 0.0;
          phi_el_old[i] = 0.0;
          phi_ed[i] = 0.0;
          phi_ed_old[i] = 0.0;
        } else if (type[i] == Li_Ctype) {
          phi_el[i] = phi_el_BC_Cat;
          phi_el_old[i] = phi_el_BC_Cat;
          phi_ed[i] = 0.0;  
          phi_ed_old[i] = 0.0;
        } else if (type[i] == Li_Atype) {
          phi_el[i] = 0.0;
          phi_el_old[i] = 0.0;
          phi_ed[i] = phi_ed_BC_anode;
          phi_ed_old[i] = phi_ed_BC_anode;
        }
        init_flag[i] = 1.0;
      }
    }
  } else {
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

void FixPotentialCC::pre_force(int vflag)
{
  updatePtrs();
  
  // Calculate hydrostatic stress for Li particles
  calculate_hydrostatic_stress();
}

/* ---------------------------------------------------------------------- */

void FixPotentialCC::post_force(int vflag)
{
  solve_cc_iteration();
}

/* ---------------------------------------------------------------------- */

void FixPotentialCC::updatePtrs()
{
  phi_el = fix_phi_el->vector_atom;
  phi_el_old = fix_phi_el_old->vector_atom;
  phi_ed = fix_phi_ed->vector_atom;
  phi_ed_old = fix_phi_ed_old->vector_atom;
  current_Li_SE = fix_current_Li_SE->vector_atom;
  hydrostatic_stress = fix_hydrostatic_stress->vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixPotentialCC::solve_cc_iteration()
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

  // Use the user-specified applied current density (from BC_potentials phi_el_BC_An)
  double i_app = phi_el_BC_An;

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

  // CC iteration for both potentials
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    
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
        double contact_area = calculate_contact_area(i, j); // m²
        
        if (contact_area > 0.0) {
          // === Electrolyte potential update (SE) ===
          if (type[i] == SE_type || type[i] == Li_Atype || type[i] == Li_Ctype) {
            if (type[i] == Li_Atype && type[j] == SE_type) {
              double conductance = sigma_el * contact_area / r_SI;
              phi_el_sum += conductance * phi_el[j];
              coeff_el_sum += conductance;
              cur_sum += i_app * contact_area;
              current_Li_SE[i] += i_app * contact_area;
            
            } else if (type[i] == Li_Atype && type[j] == Li_Atype) {
              double conductance = sigma_el * contact_area / r_SI;
              phi_el_sum += conductance * phi_el[j];
              coeff_el_sum += conductance;

            } else if (type[i] == SE_type && type[j] == SE_type) {
              double conductance = sigma_el * contact_area / r_SI;
              phi_el_sum += conductance * phi_el[j];
              coeff_el_sum += conductance;
            
            } else if (type[i] == SE_type && type[j] == Li_Atype) {
              double conductance = sigma_el * contact_area / r_SI;
              phi_el_sum += conductance * phi_el[j];
              coeff_el_sum += conductance;

            } else if (type[i] == SE_type && type[j] == Li_Ctype) {
              double conductance = sigma_el * contact_area / r_SI;
              phi_el_sum += conductance * phi_el[j];
              coeff_el_sum += conductance;
            
            } else if (type[i] == Li_Ctype && type[j] == SE_type) {
              double conductance = sigma_el * contact_area / r_SI;
              phi_el_sum += conductance * phi_el[j];
              coeff_el_sum += conductance;
              cur_sum += -1.0 * i_app * contact_area;
              current_Li_SE[i] += -1.0 * i_app * contact_area;

            } else if (type[i] == Li_Ctype && type[j] == Li_Ctype) {
              double conductance = sigma_el * contact_area / r_SI;
              phi_el_sum += conductance * phi_el[j];
              coeff_el_sum += conductance;
            }
          }
        }
      }
    }
    
    // Update electrolyte potential
    if ((type[i] == SE_type || type[i] == Li_Atype || type[i] == Li_Ctype)) {
      double phi_el_new = (phi_el_sum + cur_sum) / coeff_el_sum;

      if (!std::isnan(phi_el_new) && !std::isinf(phi_el_new)) {
        phi_el[i] = phi_el_old[i] + omega * (phi_el_new - phi_el_old[i]);
      } else {
        phi_el[i] = phi_el_old[i];
      }
    }
  }
  
  // Forward communication for both potentials
  fix_phi_el->do_forward_comm();
  fix_phi_ed->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

void FixPotentialCC::apply_boundary_conditions()
{
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (type[i] == Li_Atype) {
        // Anode boundary - potentials can be set via BC_potentials and phi_ed_anode
      }
      else if (type[i] == Li_Ctype) {
        // Cathode boundary - potentials can be set via BC_potentials
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPotentialCC::calculate_hydrostatic_stress()
{
  double **f = atom->f;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && type[i] == Li_Ctype) {
      double force_conversion = 1.0e-9;  // LAMMPS force micro nN to N SI
      double fmag = sqrt(f[i][0]*f[i][0] + f[i][1]*f[i][1] + f[i][2]*f[i][2]) * force_conversion;
      
      double radius_SI = radius[i] * 1.0e-6;
      double surface_area = 4.0 * M_PI * radius_SI * radius_SI;
      
      if (surface_area > 0.0) {
        hydrostatic_stress[i] = 0.0; // Can be enabled: fmag / surface_area
      } else {
        hydrostatic_stress[i] = 0.0;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixPotentialCC::calculate_contact_area(int i, int j)
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

double FixPotentialCC::calculate_current_Li_SE(int i_Li, int j_SE, double phi_ed_val, double phi_el_val, double sigma_m)
{
  // Use user-specified equilibrium potential and exchange current density
  double eta = phi_ed_val - phi_el_val - U_eq;

  double RT = R * T;
  
  double arg1 = (1.0 - alpha_a) * F * eta / RT - sigma_m * 9e-6 / RT;
  double arg2 = -alpha_c * F * eta / RT - sigma_m * 9e-6 / RT;
  
  double exp_term1 = exp(arg1);
  double exp_term2 = exp(arg2);
  
  // Butler-Volmer current density
  double i_pq = i_0 * (exp_term1 - exp_term2);
  
  return i_pq;
}

/* ---------------------------------------------------------------------- */

double FixPotentialCC::check_convergence()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  double local_error = 0.0;
  int count = 0;
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (type[i] == SE_type) {
        double diff_el = fabs(phi_el[i] - phi_el_old[i]);
        if (diff_el > local_error) local_error = diff_el;
      }
      count++;
    }
  }
  
  double global_error;
  MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_MAX, world);
  
  return global_error;
}

/* ---------------------------------------------------------------------- */

double FixPotentialCC::compute_scalar()
{
  return convergence_error;
}

/* ---------------------------------------------------------------------- */

double FixPotentialCC::compute_vector(int n)
{
  if (n == 0) return convergence_error;
  return 0.0;
}