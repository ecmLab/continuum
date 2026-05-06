/* ----------------------------------------------------------------------
    LIGGGHTS® DEM simulation engine

    fix_battery_eis.cpp - Refactored for full input-script configurability.

    NEW keyword arguments parsed from the fix command line
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      protocol       int     0 = Alternating Current, 1 = CCD stepping, 2 = Constant Currrent
      period         int     Half-cycle length in timesteps
      increase_step  double  CCD current increment per stage (A/m²)
      area           double  Cross-sectional area for I = i·A (m²)
      debug_freq     int     Timestep interval for debug output
      i_0            double  Exchange current density for BV (A/m²)

    Example fix lines  (see bottom of file for full commented blocks)
    ~~~~~~~~~~~~~~~~~
    # Constant-Current mode
    fix bat all battery/eis omega 1.0 tolerance 1e-1 max_iter 1 &
        BC_types 2 3 4 5 BC_potentials 0.0 10.0 &
        conductivity 0.05 0.1 SE_type 1 &
        protocol 0 area 4e-10 debug_freq 1440 i_0 0.01

    # CCD stepping mode
    fix bat all battery/eis omega 1.0 tolerance 1e-1 max_iter 1 &
        BC_types 2 3 4 5 BC_potentials 0.0 10.0 &
        conductivity 0.05 0.1 SE_type 1 &
        protocol 1 period 1800000 increase_step 2.5 &
        area 4e-10 debug_freq 720 i_0 0.01
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
#include "comm.h"
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL 1e-15
#define MAX_EXP_ARG 20.0

/* ======================================================================
   post_create  –  register all per-atom properties
   ====================================================================== */

void FixBatteryEIS::post_create()
{
  // Helper lambda-style macro to cut boilerplate
  #define REGISTER_PROP(handle, name_str)                                     \
  handle = static_cast<FixPropertyAtom*>(                                     \
      modify->find_fix_property(name_str,"property/atom","scalar",            \
                                0,0,style,false));                            \
  if (!handle) {                                                              \
    const char* fixarg[10];                                                   \
    fixarg[0] = name_str; fixarg[1] = "all";                                 \
    fixarg[2] = "property/atom"; fixarg[3] = name_str;                        \
    fixarg[4] = "scalar"; fixarg[5] = "no";                                   \
    fixarg[6] = "yes";    fixarg[7] = "no";                                   \
    fixarg[8] = "0.0";                                                        \
    handle = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);\
  }

  REGISTER_PROP(fix_phi_el,              "electrolytePotential")
  REGISTER_PROP(fix_phi_el_old,          "electrolytePotentialOld")
  REGISTER_PROP(fix_phi_ed,              "electronicPotential")
  REGISTER_PROP(fix_phi_ed_old,          "electronicPotentialOld")
  REGISTER_PROP(fix_current_Li_SE,       "currentLiSE")
  REGISTER_PROP(fix_hydrostatic_stress,  "hydrostaticStress")
  REGISTER_PROP(fix_init_flag,           "batteryInitFlag")
  REGISTER_PROP(fix_contact_area_anode,  "contactAreaAnode")
  REGISTER_PROP(fix_contact_area_cathode,"contactAreaCathode")

  #undef REGISTER_PROP
}

/* ======================================================================
   Constructor  –  parse ALL keyword/value pairs from the fix command
   ====================================================================== */

FixBatteryEIS::FixBatteryEIS(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  // --- SOR defaults ---
  omega(1.9),
  tolerance(1e-6),
  max_iterations(10000),
  current_iteration(0),
  convergence_error(0.0),
  // --- Physical constants ---
  R(8.31),
  T(303.0),
  F(96485.0),
  alpha_a(0.5),
  alpha_c(0.5),
  // --- Conductivities ---
  sigma_el(0.05),
  sigma_Li(1.0e7),
  // --- Type IDs ---
  SE_type(1),
  Li_Ctype(2),
  Li_Atype(3),
  CC_Ctype(4),
  CC_Atype(5),
  // --- BC values ---
  phi_ref(0.0),
  cur_app(1.0),
  // --- NEW configurable parameters (sensible defaults) ---
  protocol(0),              // default: Constant Current
  period(1800000),          // default half-cycle length (timesteps)
  increase_step(0.0),       // default: no stepping
  cross_area(4.0e-10),      // default: 20 µm × 20 µm
  debug_freq(720),          // default: print every 720 steps
  i_0(0.01),                // default exchange current density (A/m²)
  // --- Calculated interface values ---
  total_current(0.0),
  global_area_anode(0.0),
  global_area_cathode(0.0),
  i_density_anode(0.0),
  i_density_cathode(0.0),
  // --- Per-atom pointers ---
  phi_el(NULL),
  phi_el_old(NULL),
  phi_ed(NULL),
  phi_ed_old(NULL),
  current_Li_SE(NULL),
  hydrostatic_stress(NULL),
  contact_area_anode(NULL),
  contact_area_cathode(NULL),
  // --- Fix handles ---
  fix_phi_el(NULL),
  fix_phi_el_old(NULL),
  fix_phi_ed(NULL),
  fix_phi_ed_old(NULL),
  fix_current_Li_SE(NULL),
  fix_hydrostatic_stress(NULL),
  fix_init_flag(NULL),
  fix_contact_area_anode(NULL),
  fix_contact_area_cathode(NULL),
  // --- Misc ---
  list(NULL),
  first_run(true)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix battery/eis command: too few arguments");

  /* ------------------------------------------------------------------
     Keyword / value parsing
     ------------------------------------------------------------------ */
  int iarg = 3;
  while (iarg < narg) {

    // ---- omega (SOR relaxation factor) ----
    if (strcmp(arg[iarg],"omega") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"fix battery/eis: 'omega' requires 1 value");
      omega = force->numeric(FLERR,arg[iarg+1]);
      if (omega <= 0.0 || omega >= 2.0)
        error->all(FLERR,"fix battery/eis: omega must be in (0, 2)");
      iarg += 2;

    // ---- tolerance ----
    } else if (strcmp(arg[iarg],"tolerance") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"fix battery/eis: 'tolerance' requires 1 value");
      tolerance = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;

    // ---- max_iter ----
    } else if (strcmp(arg[iarg],"max_iter") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"fix battery/eis: 'max_iter' requires 1 value");
      max_iterations = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;

    // ---- temperature ----
    } else if (strcmp(arg[iarg],"temperature") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"fix battery/eis: 'temperature' requires 1 value");
      T = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;

    // ---- BC_types (4 atom-type IDs) ----
    } else if (strcmp(arg[iarg],"BC_types") == 0) {
      if (iarg+5 > narg)
        error->all(FLERR,"fix battery/eis: 'BC_types' requires 4 values "
                         "(Li_Ctype Li_Atype CC_Ctype CC_Atype)");
      Li_Ctype = force->inumeric(FLERR,arg[iarg+1]);
      Li_Atype = force->inumeric(FLERR,arg[iarg+2]);
      CC_Ctype = force->inumeric(FLERR,arg[iarg+3]);
      CC_Atype = force->inumeric(FLERR,arg[iarg+4]);
      iarg += 5;

    // ---- BC_potentials (phi_ref, cur_app) ----
    } else if (strcmp(arg[iarg],"BC_potentials") == 0) {
      if (iarg+3 > narg)
        error->all(FLERR,"fix battery/eis: 'BC_potentials' requires 2 values "
                         "(phi_ref cur_app)");
      phi_ref = force->numeric(FLERR,arg[iarg+1]);
      cur_app = force->numeric(FLERR,arg[iarg+2]);
      iarg += 3;

    // ---- conductivity (sigma_SE, sigma_Li) ----
    } else if (strcmp(arg[iarg],"conductivity") == 0) {
      if (iarg+3 > narg)
        error->all(FLERR,"fix battery/eis: 'conductivity' requires 2 values "
                         "(sigma_SE sigma_Li)");
      sigma_el = force->numeric(FLERR,arg[iarg+1]);
      sigma_Li = force->numeric(FLERR,arg[iarg+2]);
      iarg += 3;

    // ---- SE_type ----
    } else if (strcmp(arg[iarg],"SE_type") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"fix battery/eis: 'SE_type' requires 1 value");
      SE_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;

    /* ==============================================================
       NEW CONFIGURABLE PARAMETERS
       ============================================================== */

    // ---- protocol (0 = CC, 1 = CCD) ----
    } else if (strcmp(arg[iarg],"protocol") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"fix battery/eis: 'protocol' requires 1 value "
                         "(0 = CC, 1 = CCD)");
      protocol = force->inumeric(FLERR,arg[iarg+1]);
      if (protocol < 0 || protocol > 2)
        error->all(FLERR,"fix battery/eis: protocol must be 0 (AC), 1 (CCD), 2 (CC)");
      iarg += 2;

    // ---- period (half-cycle length in timesteps, used by CCD) ----
    } else if (strcmp(arg[iarg],"period") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"fix battery/eis: 'period' requires 1 value "
                         "(half-cycle length in timesteps)");
      period = force->inumeric(FLERR,arg[iarg+1]);
      if (period <= 0)
        error->all(FLERR,"fix battery/eis: period must be > 0");
      iarg += 2;

    // ---- increase_step (current increment per CCD stage, A/m²) ----
    } else if (strcmp(arg[iarg],"increase_step") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"fix battery/eis: 'increase_step' requires 1 value (A/m²)");
      increase_step = force->numeric(FLERR,arg[iarg+1]);
      if (increase_step < 0.0)
        error->all(FLERR,"fix battery/eis: increase_step must be >= 0");
      iarg += 2;

    // ---- area (cross-sectional area in m²) ----
    } else if (strcmp(arg[iarg],"area") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"fix battery/eis: 'area' requires 1 value (m²)");
      cross_area = force->numeric(FLERR,arg[iarg+1]);
      if (cross_area <= 0.0)
        error->all(FLERR,"fix battery/eis: area must be > 0");
      iarg += 2;

    // ---- debug_freq (timestep interval for screen output) ----
    } else if (strcmp(arg[iarg],"debug_freq") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"fix battery/eis: 'debug_freq' requires 1 value");
      debug_freq = force->inumeric(FLERR,arg[iarg+1]);
      if (debug_freq <= 0)
        error->all(FLERR,"fix battery/eis: debug_freq must be > 0");
      iarg += 2;

    // ---- i_0 (exchange current density for Butler-Volmer, A/m²) ----
    } else if (strcmp(arg[iarg],"i_0") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"fix battery/eis: 'i_0' requires 1 value (A/m²)");
      i_0 = force->numeric(FLERR,arg[iarg+1]);
      if (i_0 <= 0.0)
        error->all(FLERR,"fix battery/eis: i_0 must be > 0");
      iarg += 2;

    } else {
      error->all(FLERR,"Illegal fix battery/eis command: "
                       "unrecognised keyword");
    }
  }

  /* ------------------------------------------------------------------
     Sanity check: CCD mode requires a nonzero increase_step
     ------------------------------------------------------------------ */
  if (protocol == 1 && increase_step <= 0.0) {
    error->all(FLERR,"fix battery/eis: CCD protocol (1) requires "
                     "increase_step > 0");
  }

  /* ------------------------------------------------------------------
     Register global output flags
     ------------------------------------------------------------------ */
  scalar_flag  = 1;
  global_freq  = 1;
  extscalar    = 0;

  vector_flag  = 1;
  size_vector  = 6;
  extvector    = 0;
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

  fix_phi_el              = static_cast<FixPropertyAtom*>(modify->find_fix_property("electrolytePotential","property/atom","scalar",0,0,style));
  fix_phi_el_old          = static_cast<FixPropertyAtom*>(modify->find_fix_property("electrolytePotentialOld","property/atom","scalar",0,0,style));
  fix_phi_ed              = static_cast<FixPropertyAtom*>(modify->find_fix_property("electronicPotential","property/atom","scalar",0,0,style));
  fix_phi_ed_old          = static_cast<FixPropertyAtom*>(modify->find_fix_property("electronicPotentialOld","property/atom","scalar",0,0,style));
  fix_current_Li_SE       = static_cast<FixPropertyAtom*>(modify->find_fix_property("currentLiSE","property/atom","scalar",0,0,style));
  fix_hydrostatic_stress  = static_cast<FixPropertyAtom*>(modify->find_fix_property("hydrostaticStress","property/atom","scalar",0,0,style));
  fix_init_flag           = static_cast<FixPropertyAtom*>(modify->find_fix_property("batteryInitFlag","property/atom","scalar",0,0,style));
  fix_contact_area_anode  = static_cast<FixPropertyAtom*>(modify->find_fix_property("contactAreaAnode","property/atom","scalar",0,0,style));
  fix_contact_area_cathode= static_cast<FixPropertyAtom*>(modify->find_fix_property("contactAreaCathode","property/atom","scalar",0,0,style));

  if (!fix_phi_el || !fix_phi_el_old || !fix_phi_ed || !fix_phi_ed_old ||
      !fix_current_Li_SE || !fix_hydrostatic_stress || !fix_init_flag ||
      !fix_contact_area_anode || !fix_contact_area_cathode)
    error->all(FLERR,"Could not find required property/atom fixes");

  if (SE_type < 1 || SE_type > atom->ntypes)
    error->all(FLERR,"Invalid SE particle type for battery/eis");
  if (Li_Ctype < 1 || Li_Ctype > atom->ntypes ||
      Li_Atype < 1 || Li_Atype > atom->ntypes)
    error->all(FLERR,"Invalid Li particle types for battery/eis");
  if (CC_Ctype < 1 || CC_Ctype > atom->ntypes ||
      CC_Atype < 1 || CC_Atype > atom->ntypes)
    error->all(FLERR,"Invalid CC particle types for battery/eis");

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix  = 1;
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
  int *type  = atom->type;
  int *mask  = atom->mask;
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
        phi_el[i]     = phi_ref;
        phi_el_old[i] = phi_ref;
        phi_ed[i]     = 0.0;
        phi_ed_old[i] = 0.0;
        init_flag[i]  = 1.0;
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
  calculate_interface_currents();

  for (int iter = 0; iter < max_iterations; iter++) {
    solve_eis_iteration();
  }
}

/* ---------------------------------------------------------------------- */

void FixBatteryEIS::updatePtrs()
{
  phi_el              = fix_phi_el->vector_atom;
  phi_el_old          = fix_phi_el_old->vector_atom;
  phi_ed              = fix_phi_ed->vector_atom;
  phi_ed_old          = fix_phi_ed_old->vector_atom;
  current_Li_SE       = fix_current_Li_SE->vector_atom;
  hydrostatic_stress  = fix_hydrostatic_stress->vector_atom;
  contact_area_anode  = fix_contact_area_anode->vector_atom;
  contact_area_cathode= fix_contact_area_cathode->vector_atom;
}

/* ======================================================================
   calculate_interface_currents
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Computes contact areas at the anode and cathode Li/SE interfaces,
   then determines the applied current density using the user-selected
   protocol:

     protocol 0 (CC):  constant i_density = sign * cur_app
     protocol 1 (CCD): i_density ramps by increase_step every two
                        half-cycles

   Total current I = i_density * cross_area   (user-supplied area).
   ====================================================================== */

void FixBatteryEIS::calculate_interface_currents()
{
  int i,j,ii,jj,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,r;

  double **x   = atom->x;
  double *radius = atom->radius;
  int *type    = atom->type;
  int *mask    = atom->mask;
  int nlocal   = atom->nlocal;

  inum       = list->inum;
  ilist      = list->ilist;
  numneigh   = list->numneigh;
  firstneigh = list->firstneigh;

  // Reset per-particle contact areas
  for (i = 0; i < nlocal; i++) {
    contact_area_anode[i]  = 0.0;
    contact_area_cathode[i] = 0.0;
  }

  double local_area_anode_SE  = 0.0;
  double local_area_cathode_SE = 0.0;

  // --- First pass: accumulate per-particle and total contact areas ---
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum  = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq  = delx*delx + dely*dely + delz*delz;
      r    = sqrt(rsq);

      if (r < (radius[i] + radius[j])) {
        double ca = calculate_contact_area(i, j);

        if (ca > 0.0) {
          // Anode interface: Li_Atype <-> SE_type
          if (type[i] == Li_Atype && type[j] == SE_type) {
            contact_area_anode[i] += ca;
            local_area_anode_SE   += ca;
          } else if (type[i] == SE_type && type[j] == Li_Atype) {
            contact_area_anode[i] += ca;
          }

          // Cathode interface: Li_Ctype <-> SE_type
          if (type[i] == Li_Ctype && type[j] == SE_type) {
            contact_area_cathode[i] += ca;
            local_area_cathode_SE   += ca;
          } else if (type[i] == SE_type && type[j] == Li_Ctype) {
            contact_area_cathode[i] += ca;
          }
        }
      }
    }
  }

  // Global reduction
  MPI_Allreduce(&local_area_anode_SE,  &global_area_anode,   1,
                MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&local_area_cathode_SE, &global_area_cathode, 1,
                MPI_DOUBLE, MPI_SUM, world);

  /* ------------------------------------------------------------------
     Determine the applied current density according to protocol
     ------------------------------------------------------------------ */
  bigint ndt = update->ntimestep;
  double i_density = 0.0;

  if (protocol == 0) {
    // --- CC: simple alternating sign every half-cycle ---
    int half_cycle = static_cast<int>(ndt) / period;
    double sign = (half_cycle % 2 == 0) ? 1.0 : -1.0;
    i_density = sign * cur_app;

  } else if (protocol == 1) {
    // --- CCD: ramp magnitude by increase_step every pair of half-cycles ---
    int half_cycle       = static_cast<int>(ndt) / period;
    double sign          = (half_cycle % 2 == 0) ? 1.0 : -1.0;
    int num_pos_steps    = half_cycle / 2;
    double current_mag   = cur_app + num_pos_steps * increase_step;
    i_density = sign * current_mag;
  
  } else if (protocol == 2) {
    // --- CC: constant current ---
    i_density = cur_app;
  }

  /* ------------------------------------------------------------------
     Total current through the cell  (uses user-supplied cross_area)
     ------------------------------------------------------------------ */
  total_current = i_density * cross_area;   // Amperes

  // Distribute over actual contact areas at each interface
  i_density_anode  = (global_area_anode  > SMALL)
                   ? total_current / global_area_anode  : 0.0;
  i_density_cathode = (global_area_cathode > SMALL)
                   ? total_current / global_area_cathode : 0.0;

  /* ------------------------------------------------------------------
     Debug output at user-defined interval
     ------------------------------------------------------------------ */
  if (debug_freq > 0 && update->ntimestep % debug_freq == 0) {
    if (comm->me == 0) {
      fprintf(screen, "Battery EIS - Step %ld (protocol %d):\n",
              (long)update->ntimestep, protocol);
      fprintf(screen, "  Applied i_density = %.4f A/m²,  "
                      "cross_area = %.6e m²\n", i_density, cross_area);
      fprintf(screen, "  Anode interface:   A = %.6e m²,  "
                      "i = %.4f A/m²\n", global_area_anode, i_density_anode);
      fprintf(screen, "  Cathode interface: A = %.6e m²,  "
                      "i = %.4f A/m²\n", global_area_cathode, i_density_cathode);
      fprintf(screen, "  Total current:     I = %.6e A\n", total_current);
    }
  }
}

/* ======================================================================
   solve_eis_iteration  –  one SOR sweep
   ====================================================================== */

void FixBatteryEIS::solve_eis_iteration()
{
  int i,j,ii,jj,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,r,r_SI;
  double phi_sum,coeff_sum,cur_source;

  double **x   = atom->x;
  double *radius = atom->radius;
  int *type    = atom->type;
  int *mask    = atom->mask;
  int nlocal   = atom->nlocal;

  inum       = list->inum;
  ilist      = list->ilist;
  numneigh   = list->numneigh;
  firstneigh = list->firstneigh;

  // Store old values for convergence check
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      phi_el_old[i] = phi_el[i];
  }

  // Reset current tracking
  for (i = 0; i < nlocal; i++)
    current_Li_SE[i] = 0.0;

  // Main SOR loop
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    phi_sum    = 0.0;
    coeff_sum  = 0.0;
    cur_source = 0.0;

    jlist = firstneigh[i];
    jnum  = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq  = delx*delx + dely*dely + delz*delz;
      r    = sqrt(rsq);
      r_SI = r * 1.0e-6;

      if (r < (radius[i] + radius[j])) {
        double ca = calculate_contact_area(i, j);

        if (ca > 0.0) {
          double sigma_eff = 0.0;

          // Effective conductivity by interface type
          if (type[i] == SE_type && type[j] == SE_type) {
            sigma_eff = sigma_el;
          }
          else if ((type[i] == Li_Atype && type[j] == Li_Atype) ||
                   (type[i] == Li_Ctype && type[j] == Li_Ctype)) {
            sigma_eff = sigma_Li;
          }
          else if ((type[i] == SE_type &&
                    (type[j] == Li_Atype || type[j] == Li_Ctype)) ||
                   ((type[i] == Li_Atype || type[i] == Li_Ctype) &&
                    type[j] == SE_type)) {
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
            double conductance = sigma_eff * ca / r_SI;
            phi_sum   += conductance * phi_el[j];
            coeff_sum += conductance;
          }

          // Current source at anode (Li_Atype -> SE)
          if (type[i] == Li_Atype && type[j] == SE_type) {
            double local_cur = i_density_anode * ca;
            cur_source       += local_cur;
            current_Li_SE[i] += local_cur;
          }
          // Current sink at cathode (Li_Ctype -> SE)
          else if (type[i] == Li_Ctype && type[j] == SE_type) {
            double local_cur = -1.0 * i_density_cathode * ca;
            cur_source       += local_cur;
            current_Li_SE[i] += local_cur;
          }
        }
      }
    }

    // SOR update
    if (coeff_sum > SMALL) {
      double phi_new = (phi_sum + cur_source) / coeff_sum;

      if (!std::isnan(phi_new) && !std::isinf(phi_new)) {
        phi_el[i] = phi_el_old[i] + omega * (phi_new - phi_el_old[i]);
      } else {
        phi_el[i] = phi_el_old[i];
      }
    }
  }

  fix_phi_el->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

void FixBatteryEIS::apply_reference_potential()
{
  int *mask  = atom->mask;
  int *type  = atom->type;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && type[i] == CC_Ctype) {
      phi_el[i]     = phi_ref;
      phi_el_old[i] = phi_ref;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixBatteryEIS::calculate_hydrostatic_stress()
{
  char *stress_id = (char *)"st";
  int icompute = modify->find_compute(stress_id);

  if (icompute < 0)
    error->all(FLERR,"FixBatteryEIS: Could not find compute with ID 'st'.");

  Compute *stress_compute = modify->compute[icompute];

  if (!(stress_compute->invoked_flag & INVOKED_PERATOM)) {
    stress_compute->compute_peratom();
    stress_compute->invoked_flag |= INVOKED_PERATOM;
  }

  double **stress = stress_compute->array_atom;
  double *radius  = atom->radius;
  int *type       = atom->type;
  int *mask       = atom->mask;
  int nlocal      = atom->nlocal;

  double stress_conversion = 1.0e-15;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit &&
        (type[i] == Li_Atype || type[i] == Li_Ctype)) {
      double r_SI = radius[i] * 1.0e-6;
      double vol  = (4.0 / 3.0) * M_PI * r_SI * r_SI * r_SI;
      double trace = (stress[i][0] + stress[i][1] + stress[i][2])
                   * stress_conversion;
      hydrostatic_stress[i] = (vol > 0.0) ? trace / (3.0 * vol) : 0.0;
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixBatteryEIS::calculate_contact_area(int i, int j)
{
  double *radius = atom->radius;
  double **x     = atom->x;

  double lc = 1.0e-6;  // length conversion µm -> m
  double dx = (x[i][0] - x[j][0]) * lc;
  double dy = (x[i][1] - x[j][1]) * lc;
  double dz = (x[i][2] - x[j][2]) * lc;
  double r  = sqrt(dx*dx + dy*dy + dz*dz);

  double ri = radius[i] * lc;
  double rj = radius[j] * lc;
  double rs = ri + rj;

  if (r >= rs) return 0.0;

  if (r < fmax(ri, rj)) {
    double rmin = fmin(ri, rj);
    return M_PI * rmin * rmin;
  }

  double delta_n = rs - r;
  double reff    = ri * rj / rs;
  if (delta_n > 0.1 * fmin(ri, rj))
    delta_n = 0.1 * fmin(ri, rj);
  return M_PI * delta_n * reff;
}

/* ======================================================================
   calculate_current_Li_SE  –  Butler-Volmer kinetics
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Now uses the configurable member variable  i_0  instead of a
   hardcoded exchange current density.
   ====================================================================== */

double FixBatteryEIS::calculate_current_Li_SE(int i_Li, int j_SE,
                                               double phi_ed_val,
                                               double phi_el_val,
                                               double sigma_m)
{
  double U_eq = 0.0;
  double eta  = phi_ed_val - phi_el_val - U_eq;
  double RT   = R * T;

  double arg1 = (1.0 - alpha_a) * F * eta / RT - sigma_m * 9e-6 / RT;
  double arg2 = -alpha_c         * F * eta / RT - sigma_m * 9e-6 / RT;

  // Clamp to prevent overflow
  if (arg1 >  MAX_EXP_ARG) arg1 =  MAX_EXP_ARG;
  if (arg1 < -MAX_EXP_ARG) arg1 = -MAX_EXP_ARG;
  if (arg2 >  MAX_EXP_ARG) arg2 =  MAX_EXP_ARG;
  if (arg2 < -MAX_EXP_ARG) arg2 = -MAX_EXP_ARG;

  return i_0 * (exp(arg1) - exp(arg2));   // uses configurable i_0
}

/* ---------------------------------------------------------------------- */

double FixBatteryEIS::check_convergence()
{
  int nlocal = atom->nlocal;
  int *mask  = atom->mask;
  int *type  = atom->type;
  double local_error = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (type[i] == SE_type || type[i] == Li_Atype || type[i] == Li_Ctype) {
        double d = fabs(phi_el[i] - phi_el_old[i]);
        if (d > local_error) local_error = d;
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
  if      (n == 0) return (double)current_iteration;
  else if (n == 1) return convergence_error;
  else if (n == 2) return global_area_anode;
  else if (n == 3) return global_area_cathode;
  else if (n == 4) return i_density_anode;
  else if (n == 5) return i_density_cathode;
  return 0.0;
}