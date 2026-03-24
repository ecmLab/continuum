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
----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine
    DCS Computing GmbH, Linz, Austria

    Contributing author and copyright for this file:
    Created by: Joseph Vazquez
    Copyright 2024-     DCS Computing GmbH, Linz

    Dual-potential SOR solver for solid-state battery interlayers.
    See fix_potential_sor.h for full documentation and input-script syntax.
------------------------------------------------------------------------- */

#include "fix_potential_sor.h"
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

#define SMALL 1.0e-30

/* ======================================================================
   Helper: register a scalar property/atom if it does not already exist.
   Returns the FixPropertyAtom pointer.
   ====================================================================== */
static FixPropertyAtom *
ensure_property(Modify *modify, const char *name, const char *caller)
{
  FixPropertyAtom *fp = static_cast<FixPropertyAtom*>(
      modify->find_fix_property(name,"property/atom","scalar",0,0,caller,false));
  if (!fp) {
    const char *fixarg[9];
    fixarg[0] = name;
    fixarg[1] = "all";
    fixarg[2] = "property/atom";
    fixarg[3] = name;
    fixarg[4] = "scalar";
    fixarg[5] = "no";    // restart
    fixarg[6] = "yes";   // communicate ghost
    fixarg[7] = "no";    // communicate rev
    fixarg[8] = "0.0";
    fp = modify->add_fix_property_atom(9, const_cast<char**>(fixarg), caller);
  }
  return fp;
}

/* ======================================================================
   post_create  –  register all required property/atom fixes
   ====================================================================== */
void FixPotentialSOR::post_create()
{
  fix_phi_el      = ensure_property(modify, "electrolytePotential",    style);
  fix_phi_el_old  = ensure_property(modify, "electrolytePotentialOld", style);
  fix_phi_ed      = ensure_property(modify, "electronicPotential",     style);
  fix_phi_ed_old  = ensure_property(modify, "electronicPotentialOld",  style);
  fix_current_el  = ensure_property(modify, "currentElectrolyte",      style);
  fix_current_ed  = ensure_property(modify, "currentElectronic",       style);
  fix_hydrostatic_stress = ensure_property(modify, "hydrostaticStress",style);
  fix_init_flag   = ensure_property(modify, "batteryInitFlag",         style);
  fix_contact_area_LM = ensure_property(modify, "contactAreaLM",       style);
  fix_contact_area_SE = ensure_property(modify, "contactAreaSE",       style);
}

/* ======================================================================
   Constructor  –  parse all keyword / value arguments
   ====================================================================== */
FixPotentialSOR::FixPotentialSOR(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  /* --- SOR --- */
  omega(1.5),
  max_iterations(50),
  current_iteration(0),
  convergence_el(0.0),
  convergence_ed(0.0),
  /* --- constants --- */
  R(8.314),
  T(303.0),
  FF(96485.0),
  /* --- types (1-based) --- */
  AM_type(1), CB_type(2), LM_type(3), SE_type(4),
  /* --- geometry --- */
  area(4.0e-10),   // 20 µm × 20 µm default
  i_app(1.0),
  /* --- BC defaults: Dirichlet at LM, Neumann at SE --- */
  bc_LM_el_mode(BC_POTENTIAL), bc_LM_el_value(0.0),
  bc_LM_ed_mode(BC_POTENTIAL), bc_LM_ed_value(0.0),
  bc_SE_el_mode(BC_CURRENT),   bc_SE_el_value(1.0),  // multiplier on i_app
  bc_SE_ed_mode(BC_CURRENT),   bc_SE_ed_value(1.0),
  /* --- current tracking --- */
  total_current(0.0),
  global_area_LM(0.0), global_area_SE(0.0),
  i_dens_LM_el(0.0), i_dens_LM_ed(0.0),
  i_dens_SE_el(0.0), i_dens_SE_ed(0.0),
  /* --- per-atom pointers (set later) --- */
  phi_el(NULL), phi_el_old(NULL),
  phi_ed(NULL), phi_ed_old(NULL),
  current_el(NULL), current_ed(NULL),
  hydrostatic_stress(NULL),
  contact_area_LM(NULL), contact_area_SE(NULL),
  /* --- fix handles --- */
  fix_phi_el(NULL), fix_phi_el_old(NULL),
  fix_phi_ed(NULL), fix_phi_ed_old(NULL),
  fix_current_el(NULL), fix_current_ed(NULL),
  fix_hydrostatic_stress(NULL),
  fix_init_flag(NULL),
  fix_contact_area_LM(NULL), fix_contact_area_SE(NULL),
  /* --- misc --- */
  list(NULL),
  first_run(true)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix potential/sor command: need at least 3 args");

  // Default conductivities [S/m]:  AM       CB       LM       SE
  sigma_el_type[0] = 0.0;    // AM: no ionic conduction
  sigma_el_type[1] = 0.0;    // CB: no ionic conduction
  sigma_el_type[2] = 0.0;    // LM: set by user if desired
  sigma_el_type[3] = 0.05;   // SE: typical argyrodite

  sigma_ed_type[0] = 1.0e5;  // AM (Ag)
  sigma_ed_type[1] = 100.0;  // CB (carbon black)
  sigma_ed_type[2] = 1.0e7;  // LM
  sigma_ed_type[3] = 0.0;    // SE: negligible electronic

  /* ---- keyword parsing ---- */
  int iarg = 3;
  while (iarg < narg) {

    // ---- omega ----
    if (strcmp(arg[iarg],"omega") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"fix potential/sor: missing value after 'omega'");
      omega = force->numeric(FLERR, arg[iarg+1]);
      if (omega <= 0.0 || omega >= 2.0)
        error->all(FLERR,"fix potential/sor: omega must be in (0, 2)");
      iarg += 2;
    }
    // ---- max_iter ----
    else if (strcmp(arg[iarg],"max_iter") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"fix potential/sor: missing value after 'max_iter'");
      max_iterations = force->inumeric(FLERR, arg[iarg+1]);
      if (max_iterations < 0)
        error->all(FLERR,"fix potential/sor: max_iter must be >= 0");
      iarg += 2;
    }
    // ---- temperature ----
    else if (strcmp(arg[iarg],"temperature") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"fix potential/sor: missing value after 'temperature'");
      T = force->numeric(FLERR, arg[iarg+1]);
      if (T <= 0.0) error->all(FLERR,"fix potential/sor: temperature must be > 0");
      iarg += 2;
    }
    // ---- AM_type ----
    else if (strcmp(arg[iarg],"AM_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"fix potential/sor: missing value after 'AM_type'");
      AM_type = force->inumeric(FLERR, arg[iarg+1]);
      iarg += 2;
    }
    // ---- CB_type ----
    else if (strcmp(arg[iarg],"CB_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"fix potential/sor: missing value after 'CB_type'");
      CB_type = force->inumeric(FLERR, arg[iarg+1]);
      iarg += 2;
    }
    // ---- LM_type ----
    else if (strcmp(arg[iarg],"LM_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"fix potential/sor: missing value after 'LM_type'");
      LM_type = force->inumeric(FLERR, arg[iarg+1]);
      iarg += 2;
    }
    // ---- SE_type ----
    else if (strcmp(arg[iarg],"SE_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"fix potential/sor: missing value after 'SE_type'");
      SE_type = force->inumeric(FLERR, arg[iarg+1]);
      iarg += 2;
    }
    // ---- area ----
    else if (strcmp(arg[iarg],"area") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"fix potential/sor: missing value after 'area'");
      area = force->numeric(FLERR, arg[iarg+1]);
      if (area <= 0.0) error->all(FLERR,"fix potential/sor: area must be > 0");
      iarg += 2;
    }
    // ---- i_app ----
    else if (strcmp(arg[iarg],"i_app") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"fix potential/sor: missing value after 'i_app'");
      i_app = force->numeric(FLERR, arg[iarg+1]);
      iarg += 2;
    }
    // ---- sigma_el  <AM> <CB> <LM> <SE> ----
    else if (strcmp(arg[iarg],"sigma_el") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"fix potential/sor: sigma_el needs 4 values (AM CB LM SE)");
      for (int k = 0; k < NTYPES_MAX; k++)
        sigma_el_type[k] = force->numeric(FLERR, arg[iarg+1+k]);
      iarg += 5;
    }
    // ---- sigma_ed  <AM> <CB> <LM> <SE> ----
    else if (strcmp(arg[iarg],"sigma_ed") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"fix potential/sor: sigma_ed needs 4 values (AM CB LM SE)");
      for (int k = 0; k < NTYPES_MAX; k++)
        sigma_ed_type[k] = force->numeric(FLERR, arg[iarg+1+k]);
      iarg += 5;
    }
    // ---- BC_LM_el  <potential|current> <value> ----
    else if (strcmp(arg[iarg],"BC_LM_el") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"fix potential/sor: BC_LM_el needs mode and value");
      if      (strcmp(arg[iarg+1],"potential") == 0) bc_LM_el_mode = BC_POTENTIAL;
      else if (strcmp(arg[iarg+1],"current")  == 0) bc_LM_el_mode = BC_CURRENT;
      else error->all(FLERR,"fix potential/sor: BC_LM_el mode must be 'potential' or 'current'");
      bc_LM_el_value = force->numeric(FLERR, arg[iarg+2]);
      iarg += 3;
    }
    // ---- BC_LM_ed ----
    else if (strcmp(arg[iarg],"BC_LM_ed") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"fix potential/sor: BC_LM_ed needs mode and value");
      if      (strcmp(arg[iarg+1],"potential") == 0) bc_LM_ed_mode = BC_POTENTIAL;
      else if (strcmp(arg[iarg+1],"current")  == 0) bc_LM_ed_mode = BC_CURRENT;
      else error->all(FLERR,"fix potential/sor: BC_LM_ed mode must be 'potential' or 'current'");
      bc_LM_ed_value = force->numeric(FLERR, arg[iarg+2]);
      iarg += 3;
    }
    // ---- BC_SE_el ----
    else if (strcmp(arg[iarg],"BC_SE_el") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"fix potential/sor: BC_SE_el needs mode and value");
      if      (strcmp(arg[iarg+1],"potential") == 0) bc_SE_el_mode = BC_POTENTIAL;
      else if (strcmp(arg[iarg+1],"current")  == 0) bc_SE_el_mode = BC_CURRENT;
      else error->all(FLERR,"fix potential/sor: BC_SE_el mode must be 'potential' or 'current'");
      bc_SE_el_value = force->numeric(FLERR, arg[iarg+2]);
      iarg += 3;
    }
    // ---- BC_SE_ed ----
    else if (strcmp(arg[iarg],"BC_SE_ed") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"fix potential/sor: BC_SE_ed needs mode and value");
      if      (strcmp(arg[iarg+1],"potential") == 0) bc_SE_ed_mode = BC_POTENTIAL;
      else if (strcmp(arg[iarg+1],"current")  == 0) bc_SE_ed_mode = BC_CURRENT;
      else error->all(FLERR,"fix potential/sor: BC_SE_ed mode must be 'potential' or 'current'");
      bc_SE_ed_value = force->numeric(FLERR, arg[iarg+2]);
      iarg += 3;
    }
    /* --- legacy keywords (backward-compatible) --- */
    else if (strcmp(arg[iarg],"BC_types") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"fix potential/sor: BC_types needs 2 ints");
      LM_type = force->inumeric(FLERR, arg[iarg+1]);
      SE_type = force->inumeric(FLERR, arg[iarg+2]);
      iarg += 3;
    }
    else if (strcmp(arg[iarg],"BC_potentials") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"fix potential/sor: BC_potentials needs 3 values");
      bc_LM_el_value = force->numeric(FLERR, arg[iarg+1]);
      bc_SE_el_value = force->numeric(FLERR, arg[iarg+2]);
      i_app          = force->numeric(FLERR, arg[iarg+3]);
      iarg += 4;
    }
    else if (strcmp(arg[iarg],"conductivity") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"fix potential/sor: conductivity needs 2 values");
      sigma_el_type[3] = force->numeric(FLERR, arg[iarg+1]); // SE ionic
      sigma_ed_type[0] = force->numeric(FLERR, arg[iarg+2]); // AM electronic
      iarg += 3;
    }
    else {
      char msg[256];
      snprintf(msg, 256, "fix potential/sor: unknown keyword '%s'", arg[iarg]);
      error->all(FLERR, msg);
    }
  }

  /* ---- Sanity: at least one Dirichlet BC per field ---- */
  if (bc_LM_el_mode == BC_CURRENT && bc_SE_el_mode == BC_CURRENT)
    error->all(FLERR,"fix potential/sor: phi_el needs at least one Dirichlet BC to anchor the solution");
  if (bc_LM_ed_mode == BC_CURRENT && bc_SE_ed_mode == BC_CURRENT)
    error->all(FLERR,"fix potential/sor: phi_ed needs at least one Dirichlet BC to anchor the solution");

  /* ---- thermo output flags ---- */
  scalar_flag = 1;
  global_freq = 1;
  extscalar   = 0;
  vector_flag = 1;
  size_vector = 8;
  extvector   = 0;
}

/* ====================================================================== */
FixPotentialSOR::~FixPotentialSOR() {}

/* ====================================================================== */
int FixPotentialSOR::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  return mask;
}

/* ======================================================================
   init  –  look up fixes, validate types, request neighbour list
   ====================================================================== */
void FixPotentialSOR::init()
{
  if (force->newton_pair == 1)
    error->all(FLERR,"Fix potential/sor requires newton pair off");

  fix_phi_el      = static_cast<FixPropertyAtom*>(modify->find_fix_property("electrolytePotential",   "property/atom","scalar",0,0,style));
  fix_phi_el_old  = static_cast<FixPropertyAtom*>(modify->find_fix_property("electrolytePotentialOld","property/atom","scalar",0,0,style));
  fix_phi_ed      = static_cast<FixPropertyAtom*>(modify->find_fix_property("electronicPotential",    "property/atom","scalar",0,0,style));
  fix_phi_ed_old  = static_cast<FixPropertyAtom*>(modify->find_fix_property("electronicPotentialOld", "property/atom","scalar",0,0,style));
  fix_current_el  = static_cast<FixPropertyAtom*>(modify->find_fix_property("currentElectrolyte",     "property/atom","scalar",0,0,style));
  fix_current_ed  = static_cast<FixPropertyAtom*>(modify->find_fix_property("currentElectronic",      "property/atom","scalar",0,0,style));
  fix_hydrostatic_stress = static_cast<FixPropertyAtom*>(modify->find_fix_property("hydrostaticStress","property/atom","scalar",0,0,style));
  fix_init_flag   = static_cast<FixPropertyAtom*>(modify->find_fix_property("batteryInitFlag",        "property/atom","scalar",0,0,style));
  fix_contact_area_LM = static_cast<FixPropertyAtom*>(modify->find_fix_property("contactAreaLM",      "property/atom","scalar",0,0,style));
  fix_contact_area_SE = static_cast<FixPropertyAtom*>(modify->find_fix_property("contactAreaSE",      "property/atom","scalar",0,0,style));

  if (!fix_phi_el || !fix_phi_el_old || !fix_phi_ed || !fix_phi_ed_old ||
      !fix_current_el || !fix_current_ed || !fix_hydrostatic_stress ||
      !fix_init_flag || !fix_contact_area_LM || !fix_contact_area_SE)
    error->all(FLERR,"fix potential/sor: could not find required property/atom fixes");

  // Validate types
  int nt = atom->ntypes;
  if (AM_type < 1 || AM_type > nt || CB_type < 1 || CB_type > nt ||
      LM_type < 1 || LM_type > nt || SE_type < 1 || SE_type > nt)
    error->all(FLERR,"fix potential/sor: AM_type/CB_type/LM_type/SE_type out of range");

  // Neighbour list: full list
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix  = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  updatePtrs();
}

/* ====================================================================== */
void FixPotentialSOR::init_list(int /*id*/, NeighList *ptr) { list = ptr; }

/* ======================================================================
   setup  –  initialise potentials (once) or carry forward restart values
   ====================================================================== */
void FixPotentialSOR::setup(int /*vflag*/)
{
  updatePtrs();

  int nlocal = atom->nlocal;
  int *type  = atom->type;
  int *mask  = atom->mask;
  double *init_flag = fix_init_flag->vector_atom;

  // Detect whether particles were already initialised (e.g. from restart)
  bool already = false;
  for (int i = 0; i < nlocal; i++)
    if ((mask[i] & groupbit) && init_flag[i] > 0.5) { already = true; break; }

  if (!already) {
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;

      // Sensible initial guess: set each particle to the Dirichlet value
      // of the nearest boundary, or zero.
      double init_el = 0.0, init_ed = 0.0;

      if (is_LM(type[i])) {
        init_el = (bc_LM_el_mode == BC_POTENTIAL) ? bc_LM_el_value : 0.0;
        init_ed = (bc_LM_ed_mode == BC_POTENTIAL) ? bc_LM_ed_value : 0.0;
      } else if (is_SE(type[i])) {
        init_el = (bc_SE_el_mode == BC_POTENTIAL) ? bc_SE_el_value : 0.0;
        init_ed = (bc_SE_ed_mode == BC_POTENTIAL) ? bc_SE_ed_value : 0.0;
      }
      // Interlayer particles start at zero (will be solved)

      phi_el[i] = init_el;  phi_el_old[i] = init_el;
      phi_ed[i] = init_ed;  phi_ed_old[i] = init_ed;

      init_flag[i] = 1.0;
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

  apply_boundary_conditions();
}

/* ====================================================================== */
void FixPotentialSOR::pre_force(int /*vflag*/)
{
  updatePtrs();
  // Request per-atom virial so that compute stress/atom is available
  update->vflag_atom = update->ntimestep;
}

/* ======================================================================
   post_force  –  main solve: hydrostatic stress, then SOR iterations
   ====================================================================== */
void FixPotentialSOR::post_force(int /*vflag*/)
{
  // 1. Hydrostatic stress (needed for Butler-Volmer extensions)
  calculate_hydrostatic_stress();
  fix_hydrostatic_stress->do_forward_comm();

  // 2. Synchronise ghost potentials before iterating
  fix_phi_el->do_forward_comm();
  fix_phi_ed->do_forward_comm();

  // 3. SOR iterations
  for (int iter = 0; iter < max_iterations; iter++) {

    // Recompute interface areas and current densities each iteration
    calculate_interface_areas();

    // Compute per-interface Neumann current densities from i_app
    //
    //   Sign convention for source terms:
    //     POSITIVE = current flowing INTO the particle being updated.
    //
    //   Default physics:
    //     phi_el – ionic current enters interlayer from SE  → +i on interlayer at SE
    //     phi_ed – electronic current enters interlayer from LM → +i on interlayer at LM
    //     (current leaving at the opposite interface is handled by the
    //      Dirichlet BC, so the user only specifies one Neumann BC.)

    total_current = i_app * area;  // [A]

    // -- Electrolyte ------------------------------------------------
    if (bc_LM_el_mode == BC_CURRENT) {
      i_dens_LM_el = (global_area_LM > SMALL)
                      ? (bc_LM_el_value * total_current) / global_area_LM
                      : 0.0;
    } else {
      i_dens_LM_el = 0.0;  // Dirichlet – no explicit source
    }
    if (bc_SE_el_mode == BC_CURRENT) {
      i_dens_SE_el = (global_area_SE > SMALL)
                      ? (bc_SE_el_value * total_current) / global_area_SE
                      : 0.0;
    } else {
      i_dens_SE_el = 0.0;
    }

    // -- Electronic -------------------------------------------------
    if (bc_LM_ed_mode == BC_CURRENT) {
      i_dens_LM_ed = (global_area_LM > SMALL)
                      ? (bc_LM_ed_value * total_current) / global_area_LM
                      : 0.0;
    } else {
      i_dens_LM_ed = 0.0;
    }
    if (bc_SE_ed_mode == BC_CURRENT) {
      i_dens_SE_ed = (global_area_SE > SMALL)
                      ? (bc_SE_ed_value * total_current) / global_area_SE
                      : 0.0;
    } else {
      i_dens_SE_ed = 0.0;
    }

    // Solve electrolyte potential
    solve_potential_iteration(
        phi_el, phi_el_old,
        sigma_el_type,
        bc_LM_el_mode, bc_LM_el_value,
        bc_SE_el_mode, bc_SE_el_value,
        i_dens_LM_el,  i_dens_SE_el,
        current_el,
        fix_phi_el);

    // Solve electronic potential
    solve_potential_iteration(
        phi_ed, phi_ed_old,
        sigma_ed_type,
        bc_LM_ed_mode, bc_LM_ed_value,
        bc_SE_ed_mode, bc_SE_ed_value,
        i_dens_LM_ed,  i_dens_SE_ed,
        current_ed,
        fix_phi_ed);

    current_iteration = iter + 1;
  }

  // Store convergence metrics
  convergence_el = check_convergence(phi_el, phi_el_old);
  convergence_ed = check_convergence(phi_ed, phi_ed_old);
}

/* ======================================================================
   updatePtrs  –  refresh bare pointers after neighbour-list rebuilds
   ====================================================================== */
void FixPotentialSOR::updatePtrs()
{
  phi_el     = fix_phi_el->vector_atom;
  phi_el_old = fix_phi_el_old->vector_atom;
  phi_ed     = fix_phi_ed->vector_atom;
  phi_ed_old = fix_phi_ed_old->vector_atom;
  current_el = fix_current_el->vector_atom;
  current_ed = fix_current_ed->vector_atom;
  hydrostatic_stress = fix_hydrostatic_stress->vector_atom;
  contact_area_LM    = fix_contact_area_LM->vector_atom;
  contact_area_SE    = fix_contact_area_SE->vector_atom;
}

/* ======================================================================
   calculate_interface_areas
   Accumulates per-particle and global contact areas at the two interfaces.
   Only counts from the interlayer side to avoid double-counting.
   ====================================================================== */
void FixPotentialSOR::calculate_interface_areas()
{
  double **x      = atom->x;
  double  *radius = atom->radius;
  int     *type   = atom->type;
  int      nlocal = atom->nlocal;

  int  inum  = list->inum;
  int *ilist = list->ilist;
  int *numneigh    = list->numneigh;
  int **firstneigh = list->firstneigh;

  // Zero per-particle areas
  for (int i = 0; i < nlocal; i++) {
    contact_area_LM[i] = 0.0;
    contact_area_SE[i] = 0.0;
  }

  double local_area_LM = 0.0;
  double local_area_SE = 0.0;

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    double xtmp = x[i][0], ytmp = x[i][1], ztmp = x[i][2];
    int *jlist = firstneigh[i];
    int  jnum  = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj] & NEIGHMASK;

      double dx = xtmp - x[j][0];
      double dy = ytmp - x[j][1];
      double dz = ztmp - x[j][2];
      double r  = sqrt(dx*dx + dy*dy + dz*dz);

      if (r >= radius[i] + radius[j]) continue;

      double ca = calculate_contact_area(i, j);
      if (ca <= 0.0) continue;

      // Count from interlayer side only (avoids double-counting in MPI sum)
      if (is_interlayer(type[i]) && is_LM(type[j])) {
        contact_area_LM[i] += ca;
        local_area_LM += ca;
      }
      else if (is_LM(type[i]) && is_interlayer(type[j])) {
        contact_area_LM[i] += ca;
        // Don't add to total – counted from the interlayer side
      }

      if (is_interlayer(type[i]) && is_SE(type[j])) {
        contact_area_SE[i] += ca;
        local_area_SE += ca;
      }
      else if (is_SE(type[i]) && is_interlayer(type[j])) {
        contact_area_SE[i] += ca;
      }
    }
  }

  MPI_Allreduce(&local_area_LM, &global_area_LM, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&local_area_SE, &global_area_SE, 1, MPI_DOUBLE, MPI_SUM, world);
}

/* ======================================================================
   solve_potential_iteration  –  one SOR sweep for a single potential field
   
   Parameters
   ----------
   phi, phi_old   : per-atom current and previous potential
   sigma_type     : conductivity array [AM,CB,LM,SE]
   bc_LM_mode/val : BC at the interlayer–LM interface
   bc_SE_mode/val : BC at the interlayer–SE interface
   i_dens_LM/SE   : Neumann current density [A/m^2] at each interface
                     (0 when the corresponding BC is Dirichlet)
   current_arr    : per-atom array to store interface currents
   fix_phi        : handle for forward communication
   ====================================================================== */
void FixPotentialSOR::solve_potential_iteration(
    double *phi, double *phi_old,
    const double *sigma_type,
    BCMode bc_LM_mode, double bc_LM_val,
    BCMode bc_SE_mode, double bc_SE_val,
    double i_dens_LM, double i_dens_SE,
    double *current_arr,
    FixPropertyAtom *fix_phi)
{
  double **x      = atom->x;
  double  *radius = atom->radius;
  int     *type   = atom->type;
  int     *mask   = atom->mask;
  int      nlocal = atom->nlocal;

  int  inum  = list->inum;
  int *ilist = list->ilist;
  int *numneigh    = list->numneigh;
  int **firstneigh = list->firstneigh;

  // Save old values for this sub-iteration
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      phi_old[i] = phi[i];

  // Zero current tracking
  for (int i = 0; i < nlocal; i++)
    current_arr[i] = 0.0;

  // --- Main SOR loop over local atoms ---
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];

    // Skip Dirichlet BC particles – their potential is fixed
    if (is_LM(type[i]) && bc_LM_mode == BC_POTENTIAL) continue;
    if (is_SE(type[i]) && bc_SE_mode == BC_POTENTIAL) continue;

    double xtmp = x[i][0], ytmp = x[i][1], ztmp = x[i][2];
    int *jlist = firstneigh[i];
    int  jnum  = numneigh[i];

    double phi_sum    = 0.0;
    double coeff_sum  = 0.0;
    double cur_source = 0.0;

    int ti = type[i];
    int idx_i = type_index(ti);  // -1 for unknown types (walls, etc.)

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj] & NEIGHMASK;

      double dx = xtmp - x[j][0];
      double dy = ytmp - x[j][1];
      double dz = ztmp - x[j][2];
      double rsq = dx*dx + dy*dy + dz*dz;
      double r   = sqrt(rsq);

      if (r >= radius[i] + radius[j]) continue;

      double ca = calculate_contact_area(i, j);
      if (ca <= 0.0) continue;

      double r_SI = r * 1.0e-6;  // µm → m
      int tj    = type[j];
      int idx_j = type_index(tj);

      // ------ Determine treatment for this pair ------

      bool pair_is_LM_iface = is_LM_interface(ti, tj);
      bool pair_is_SE_iface = is_SE_interface(ti, tj);

      if (pair_is_LM_iface) {
        /*  Interlayer–LM interface  */
        if (bc_LM_mode == BC_POTENTIAL) {
          // Dirichlet: include conductance to the fixed-potential particle.
          // The LM particle's potential is clamped to bc_LM_val.
          double sig_i = (idx_i >= 0) ? sigma_type[idx_i] : 0.0;
          double sig_j = (idx_j >= 0) ? sigma_type[idx_j] : 0.0;
          double sig_eff = harmonic_mean(sig_i, sig_j);
          if (sig_eff > SMALL && r_SI > SMALL) {
            double G = sig_eff * ca / r_SI;
            // j is the LM particle → use fixed BC value
            if (is_LM(tj)) {
              phi_sum   += G * bc_LM_val;
            } else {
              // i is LM (shouldn't reach here because we skip Dirichlet particles above)
              phi_sum   += G * phi[j];
            }
            coeff_sum += G;
          }
        } else {
          // Neumann: inject current, no bulk conductance across interface
          if (is_interlayer(ti) && is_LM(tj)) {
            double local_I = i_dens_LM * ca;   // [A]
            cur_source      += local_I;
            current_arr[i]  += local_I;
          }
          else if (is_LM(ti) && is_interlayer(tj)) {
            double local_I = -i_dens_LM * ca;  // opposite sign on LM side
            cur_source      += local_I;
            current_arr[i]  += local_I;
          }
        }
      }
      else if (pair_is_SE_iface) {
        /*  Interlayer–SE interface  */
        if (bc_SE_mode == BC_POTENTIAL) {
          double sig_i = (idx_i >= 0) ? sigma_type[idx_i] : 0.0;
          double sig_j = (idx_j >= 0) ? sigma_type[idx_j] : 0.0;
          double sig_eff = harmonic_mean(sig_i, sig_j);
          if (sig_eff > SMALL && r_SI > SMALL) {
            double G = sig_eff * ca / r_SI;
            if (is_SE(tj)) {
              phi_sum += G * bc_SE_val;
            } else {
              phi_sum += G * phi[j];
            }
            coeff_sum += G;
          }
        } else {
          if (is_interlayer(ti) && is_SE(tj)) {
            double local_I = i_dens_SE * ca;
            cur_source      += local_I;
            current_arr[i]  += local_I;
          }
          else if (is_SE(ti) && is_interlayer(tj)) {
            double local_I = -i_dens_SE * ca;
            cur_source      += local_I;
            current_arr[i]  += local_I;
          }
        }
      }
      else {
        /*  Bulk pair (same phase or within interlayer)  */
        if (idx_i < 0 || idx_j < 0) continue;  // skip wall particles, etc.
        double sig_eff = harmonic_mean(sigma_type[idx_i], sigma_type[idx_j]);
        if (sig_eff > SMALL && r_SI > SMALL) {
          double G = sig_eff * ca / r_SI;
          phi_sum   += G * phi[j];
          coeff_sum += G;
        }
      }
    } // end neighbour loop

    // SOR update
    if (coeff_sum > SMALL) {
      double phi_new = (phi_sum + cur_source) / coeff_sum;
      if (std::isfinite(phi_new)) {
        phi[i] = phi_old[i] + omega * (phi_new - phi_old[i]);
      }
      // else: keep old value (numerical safety)
    }
  } // end atom loop

  // Communicate updated potentials to ghost atoms
  fix_phi->do_forward_comm();
}

/* ======================================================================
   apply_boundary_conditions  –  clamp Dirichlet particles
   ====================================================================== */
void FixPotentialSOR::apply_boundary_conditions()
{
  int *mask  = atom->mask;
  int *type  = atom->type;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    if (is_LM(type[i])) {
      if (bc_LM_el_mode == BC_POTENTIAL) {
        phi_el[i] = bc_LM_el_value;
        phi_el_old[i] = bc_LM_el_value;
      }
      if (bc_LM_ed_mode == BC_POTENTIAL) {
        phi_ed[i] = bc_LM_ed_value;
        phi_ed_old[i] = bc_LM_ed_value;
      }
    }
    else if (is_SE(type[i])) {
      if (bc_SE_el_mode == BC_POTENTIAL) {
        phi_el[i] = bc_SE_el_value;
        phi_el_old[i] = bc_SE_el_value;
      }
      if (bc_SE_ed_mode == BC_POTENTIAL) {
        phi_ed[i] = bc_SE_ed_value;
        phi_ed_old[i] = bc_SE_ed_value;
      }
    }
  }

  fix_phi_el->do_forward_comm();
  fix_phi_ed->do_forward_comm();
}

/* ======================================================================
   calculate_hydrostatic_stress  –  from compute stress/atom (ID "st")
   ====================================================================== */
void FixPotentialSOR::calculate_hydrostatic_stress()
{
  int icompute = modify->find_compute((char*)"st");
  if (icompute < 0)
    error->all(FLERR,"FixPotentialSOR: could not find compute 'st'. "
               "Add 'compute st all stress/atom' before this fix.");

  Compute *sc = modify->compute[icompute];
  if (!(sc->invoked_flag & INVOKED_PERATOM)) {
    sc->compute_peratom();
    sc->invoked_flag |= INVOKED_PERATOM;
  }

  double **stress = sc->array_atom;
  double *radius  = atom->radius;
  int    *mask    = atom->mask;
  int     nlocal  = atom->nlocal;

  // 1 pg·µm²/µs² = 1e-15 J
  const double conv = 1.0e-15;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    double r_SI = radius[i] * 1.0e-6;
    double vol  = (4.0/3.0) * M_PI * r_SI * r_SI * r_SI;
    if (vol > 0.0) {
      double trace = (stress[i][0] + stress[i][1] + stress[i][2]) * conv;
      hydrostatic_stress[i] = trace / (3.0 * vol);
    } else {
      hydrostatic_stress[i] = 0.0;
    }
  }
}

/* ======================================================================
   calculate_contact_area  –  Hertz overlap area in SI [m^2]
   ====================================================================== */
double FixPotentialSOR::calculate_contact_area(int i, int j)
{
  double **x      = atom->x;
  double  *radius = atom->radius;
  const double L  = 1.0e-6;   // µm → m

  double dx = (x[i][0] - x[j][0]) * L;
  double dy = (x[i][1] - x[j][1]) * L;
  double dz = (x[i][2] - x[j][2]) * L;
  double r  = sqrt(dx*dx + dy*dy + dz*dz);

  double ri = radius[i] * L;
  double rj = radius[j] * L;
  double rsum = ri + rj;

  if (r >= rsum) return 0.0;

  if (r < fmax(ri, rj)) {
    // One sphere inside the other
    double rmin = fmin(ri, rj);
    return M_PI * rmin * rmin;
  }

  // Hertz contact area
  double delta = rsum - r;
  double reff  = ri * rj / rsum;
  // Cap overlap at 10% of smallest radius to avoid unphysical areas
  double dmax  = 0.1 * fmin(ri, rj);
  if (delta > dmax) delta = dmax;
  return M_PI * delta * reff;
}

/* ======================================================================
   check_convergence  –  max |phi - phi_old| across all local atoms
   ====================================================================== */
double FixPotentialSOR::check_convergence(double *phi, double *phi_old)
{
  int nlocal = atom->nlocal;
  int *mask  = atom->mask;
  double local_err = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double d = fabs(phi[i] - phi_old[i]);
      if (d > local_err) local_err = d;
    }
  }
  double global_err;
  MPI_Allreduce(&local_err, &global_err, 1, MPI_DOUBLE, MPI_MAX, world);
  return global_err;
}

/* ======================================================================
   thermo output
   ====================================================================== */
double FixPotentialSOR::compute_scalar()
{
  // Return the larger of the two convergence errors
  return fmax(convergence_el, convergence_ed);
}

double FixPotentialSOR::compute_vector(int n)
{
  switch (n) {
    case 0: return (double)current_iteration;
    case 1: return convergence_el;
    case 2: return convergence_ed;
    case 3: return global_area_LM;
    case 4: return global_area_SE;
    case 5: return i_dens_LM_el + i_dens_LM_ed;  // total current density at LM
    case 6: return i_dens_SE_el + i_dens_SE_ed;  // total current density at SE
    case 7: return total_current;
    default: return 0.0;
  }
}