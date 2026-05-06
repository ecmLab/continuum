/* ----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine
    Contributing author: Joseph Vazquez Mercado, RIT 2025
    Copyright 2024-     DCS Computing GmbH, Linz
    
    Implicit (Backward Euler) Lithium Diffusion on Particle Contact Network
    -----------------------------------------------------------------------
    Replaces explicit Euler + scaling factor with an unconditionally stable
    implicit solver using Jacobi iteration.
    
    - Same-type pairs (AM-AM, CB-CB): IMPLICIT concentration-driven (Fick)
    - Cross-type pairs (AM-CB, AM-LM, CB-AM, CB-LM): SEMI-IMPLICIT
      chemo-mechanical driving force:
        Φ_eff = U_OCV - (σ_h * Ω_Li) / F
      where σ_h is hydrostatic stress and Ω_Li is the partial molar volume
    - LM particles: Dirichlet BC (constant concentration)
    
    The implicit system is diagonally dominant → Jacobi convergence guaranteed.
    
    Diffusion coefficients:
      - AM: constant (user-specified, default 1.0e-10 m²/s)
      - LM: constant (user-specified, default 1.0e-14 m²/s)
      - CB: SOC-dependent via log-space rational fit from GITT data
            log10(D_CB) = P(x)/Q(x),  D_CB = 10^(P/Q)
      - Pair coefficient: min(D_i, D_j) — the limiting species
    
    Current coupling:
      The electrolyte current (current_el) computed by fix potential/sor
      ("currentElectrolyte" property/atom) is consumed here as an explicit
      source term in the cross-flux accumulation.
      
    Usage in input script:
      fix ID group lithium_diffusion_implicit AM_type 1 CB_type 2 LM_type 3 \
          diffusion_dt 0.5 update_every 100 max_iter 200 tolerance 1e-10 \
          c_li_max_LM 77101.002 c_li_max_AM 85155.85 c_li_max_CB 27776.39 \
          V_exp_max_AM 10.291 \
          D_AM 1.0e-10 D_LM 1.0e-14 \
          Omega_Li_AM 1.0602e-05 Omega_Li_CB 4.1418e-06 Omega_Li_LM 12.97e-06
------------------------------------------------------------------------- */

#include "fix_lithium_diffusion_implicit.h"
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
#include "memory.h"
#include "comm.h"
#include <cmath>
#include <cstring>
#include <algorithm>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixLithiumDiffusionImplicit::FixLithiumDiffusionImplicit(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  R(8.314),
  T(303.0),
  F(96485.0),
  c_li_max_LM(77101.002),
  c_li_max_AM(85155.85),
  c_li_max_CB(23332.2),
  V_exp_max_AM(1.01),
  V_exp_max_CB(1.13),
  V_exp_max_LM(0.0),
  Omega_Li_AM(1.0602e-05),
  Omega_Li_CB(4.1418e-06),
  Omega_Li_LM(12.97e-06),
  initial_lithium_content(0.0),
  target_lithium_content(1.0),
  max_lithium_content(1.0),
  D_AM(1.0e-10),
  D_LM(1.0e-14),
  D_CB_floor(1.0e-20),
  diffusion_dt(0.1),
  update_every(100),
  max_iter(200),
  tol(1.0e-10),
  num_substeps(1),
  nmax_alloc(0),
  c_old(NULL),
  c_new(NULL),
  cross_flux_arr(NULL),
  lithium_content(NULL),
  lithium_concentration(NULL),
  equilibrium_potential(NULL),
  hydrostatic_stress(NULL),
  current_el(NULL),
  diffusion_coefficient(NULL),
  lithium_flux(NULL),
  li_mols(NULL),
  initial_volume(NULL),
  fix_initial_volume(NULL),
  fix_lithium_content(NULL),
  fix_lithium_concentration(NULL),
  fix_equilibrium_potential(NULL),
  fix_hydrostatic_stress(NULL),
  fix_current_el(NULL),
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
    error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"AM_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      AM_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"CB_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      CB_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"LM_type") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      LM_type = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"temperature") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      T = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"diffusion_dt") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      diffusion_dt = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"update_every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      update_every = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"max_iter") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      max_iter = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"tolerance") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      tol = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"num_substeps") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      num_substeps = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"c_li_max_LM") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      c_li_max_LM = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"c_li_max_AM") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      c_li_max_AM = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"c_li_max_CB") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      c_li_max_CB = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"D_AM") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      D_AM = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"D_LM") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      D_LM = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"D_CB_floor") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      D_CB_floor = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"V_exp_max_AM") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      V_exp_max_AM = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"V_exp_max_CB") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      V_exp_max_CB = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"V_exp_max_LM") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      V_exp_max_LM = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"Omega_Li_AM") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      Omega_Li_AM = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"Omega_Li_CB") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      Omega_Li_CB = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"Omega_Li_LM") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
      Omega_Li_LM = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix lithium_diffusion_implicit command");
  }

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;

  vector_flag = 1;
  size_vector = 3;
  extvector = 0;

  last_iter_count = 0;
  last_residual = 0.0;
  total_diffusion_time = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusionImplicit::post_create()
{
  fix_diffusion_coefficient = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("diffusionCoefficient","property/atom","scalar",0,0,style,false));
  if(!fix_diffusion_coefficient) {
    const char* fixarg[10];
    fixarg[0]="diffusionCoefficient"; fixarg[1]="all"; fixarg[2]="property/atom";
    fixarg[3]="diffusionCoefficient"; fixarg[4]="scalar"; fixarg[5]="no";
    fixarg[6]="yes"; fixarg[7]="no"; fixarg[8]="0.0";
    fix_diffusion_coefficient = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  fix_lithium_flux = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("lithiumFlux","property/atom","scalar",0,0,style,false));
  if(!fix_lithium_flux) {
    const char* fixarg[10];
    fixarg[0]="lithiumFlux"; fixarg[1]="all"; fixarg[2]="property/atom";
    fixarg[3]="lithiumFlux"; fixarg[4]="scalar"; fixarg[5]="no";
    fixarg[6]="yes"; fixarg[7]="no"; fixarg[8]="0.0";
    fix_lithium_flux = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  fix_li_mols = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("lithiumMols","property/atom","scalar",0,0,style,false));
  if(!fix_li_mols) {
    const char* fixarg[10];
    fixarg[0]="lithiumMols"; fixarg[1]="all"; fixarg[2]="property/atom";
    fixarg[3]="lithiumMols"; fixarg[4]="scalar"; fixarg[5]="no";
    fixarg[6]="yes"; fixarg[7]="no"; fixarg[8]="0.0";
    fix_li_mols = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  fix_lithium_concentration = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("lithiumConcentration","property/atom","scalar",0,0,style,false));
  if(!fix_lithium_concentration) {
    const char* fixarg[10];
    fixarg[0]="lithiumConcentration"; fixarg[1]="all"; fixarg[2]="property/atom";
    fixarg[3]="lithiumConcentration"; fixarg[4]="scalar"; fixarg[5]="no";
    fixarg[6]="yes"; fixarg[7]="no"; fixarg[8]="0.0";
    fix_lithium_concentration = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  fix_initial_volume = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("initialVolume","property/atom","scalar",0,0,style,false));
  if(!fix_initial_volume) {
    const char* fixarg[10];
    fixarg[0]="initialVolume"; fixarg[1]="all"; fixarg[2]="property/atom";
    fixarg[3]="initialVolume"; fixarg[4]="scalar"; fixarg[5]="no";
    fixarg[6]="yes"; fixarg[7]="no"; fixarg[8]="0.0";
    fix_initial_volume = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }
}

/* ---------------------------------------------------------------------- */

FixLithiumDiffusionImplicit::~FixLithiumDiffusionImplicit()
{
  memory->destroy(c_old);
  memory->destroy(c_new);
  memory->destroy(cross_flux_arr);
}

/* ---------------------------------------------------------------------- */

int FixLithiumDiffusionImplicit::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusionImplicit::init()
{
  fix_lithium_content = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("lithiumContent","property/atom","scalar",0,0,style));
  if(!fix_lithium_content)
    error->all(FLERR,"Fix lithium_diffusion_implicit requires lithiumContent property");

  fix_lithium_concentration = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("lithiumConcentration","property/atom","scalar",0,0,style));
  if(!fix_lithium_concentration)
    error->all(FLERR,"Fix lithium_diffusion_implicit requires lithiumConcentration property");

  fix_initial_volume = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("initialVolume","property/atom","scalar",0,0,style));
  if(!fix_initial_volume)
    error->all(FLERR,"Could not find initialVolume property/atom fix");

  fix_equilibrium_potential = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("equilibriumPotential","property/atom","scalar",0,0,style));
  if(!fix_equilibrium_potential)
    error->all(FLERR,"Fix lithium_diffusion_implicit requires equilibriumPotential property");

  fix_hydrostatic_stress = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("hydrostaticStress","property/atom","scalar",0,0,style));
  if(!fix_hydrostatic_stress)
    error->all(FLERR,"Fix lithium_diffusion_implicit requires hydrostaticStress property "
               "(created by fix potential/sor)");

  fix_current_el = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("currentElectrolyte","property/atom","scalar",0,0,style));
  if(!fix_current_el)
    error->all(FLERR,"Fix lithium_diffusion_implicit requires currentElectrolyte property "
               "(created by fix potential/sor)");

  fix_diffusion_coefficient = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("diffusionCoefficient","property/atom","scalar",0,0,style));
  fix_lithium_flux = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("lithiumFlux","property/atom","scalar",0,0,style));
  fix_li_mols = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("lithiumMols","property/atom","scalar",0,0,style));

  if(!fix_diffusion_coefficient || !fix_lithium_flux || !fix_li_mols)
    error->all(FLERR,"Could not find required property/atom fixes");

  // Find lithium content manager for initial/target/max values
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

  // Request full neighbor list
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusionImplicit::setup(int vflag)
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

  // Initialize per-atom diffusion coefficients
  update_per_atom_D();

  total_diffusion_time = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusionImplicit::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusionImplicit::end_of_step()
{
  if (update->ntimestep % update_every != 0) return;

  updatePtrs();
  ensure_arrays_sized();

  // Synchronize ghost atoms before the solve
  fix_lithium_content->do_forward_comm();
  fix_lithium_concentration->do_forward_comm();
  fix_equilibrium_potential->do_forward_comm();
  fix_hydrostatic_stress->do_forward_comm();
  fix_current_el->do_forward_comm();

  // Update per-atom D from current lithium content (before sub-stepping)
  update_per_atom_D();
  fix_diffusion_coefficient->do_forward_comm();

  double sub_dt = diffusion_dt / num_substeps;

  for (int ss = 0; ss < num_substeps; ss++) {
    implicit_diffusion_step(sub_dt);
  }

  total_diffusion_time += diffusion_dt;

  // Forward communicate all updated quantities
  fix_lithium_content->do_forward_comm();
  fix_lithium_concentration->do_forward_comm();
  fix_li_mols->do_forward_comm();
}

/* ---------------------------------------------------------------------- */
/*  Update the per-atom diffusionCoefficient property from current
 *  lithium content.  AM and LM get constants; CB gets the GITT fit.
 * ---------------------------------------------------------------------- */

void FixLithiumDiffusionImplicit::update_per_atom_D()
{
  int nlocal = atom->nlocal;
  int *type  = atom->type;
  int *mask  = atom->mask;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    diffusion_coefficient[i] = compute_D_particle(i);
  }
}

/* ---------------------------------------------------------------------- */
/*  Compute the diffusion coefficient for a single particle.
 *    AM  → constant D_AM
 *    LM  → constant D_LM
 *    CB  → SOC-dependent from GITT rational fit
 * ---------------------------------------------------------------------- */

double FixLithiumDiffusionImplicit::compute_D_particle(int i)
{
  int t = atom->type[i];

  if (t == AM_type) return D_AM;
  if (t == LM_type) return D_LM;

  if (t == CB_type) {
    double x_norm = lithium_content[i] / max_lithium_content;
    return compute_D_CB(x_norm);
  }

  return 0.0;
}

/* ---------------------------------------------------------------------- */
/*  Carbon black diffusion coefficient from GITT log-space rational fit.
 *
 *  METHOD A (R²_log = 0.987):
 *    log10(D) = P(x) / Q(x)
 *    P(x) = 9th-order polynomial in x_norm
 *    Q(x) = 8th-order polynomial in x_norm
 *    D    = 10^( P/Q )
 *
 *  x_norm is clamped to [1e-6, 1-1e-6] to avoid singularities.
 *  Result is clamped to [D_CB_floor, 1e-8] for physical sanity.
 * ---------------------------------------------------------------------- */

double FixLithiumDiffusionImplicit::compute_D_CB(double x_norm)
{
  // Clamp to avoid singularities at boundaries
  if (x_norm < 1.0e-6)  x_norm = 1.0e-6;
  if (x_norm > 1.0 - 1.0e-6) x_norm = 1.0 - 1.0e-6;

  double x2 = x_norm * x_norm;
  double x3 = x2 * x_norm;
  double x4 = x3 * x_norm;
  double x5 = x4 * x_norm;
  double x6 = x5 * x_norm;
  double x7 = x6 * x_norm;
  double x8 = x7 * x_norm;
  double x9 = x8 * x_norm;

  // Numerator P(x) — 9th order
  double P = -18.1316386996 * x9
           +  61.5176399234 * x8
           +   3.7451520749 * x7
           -  50.4131104634 * x6
           -  16.6597657289 * x5
           +  49.8491805574 * x4
           -  25.2979337653 * x3
           +   5.6168931983 * x2
           -   0.5769340992 * x_norm
           +   0.0224398257;

  // Denominator Q(x) — 8th order
  double Q =  1.0            * x8
           -  9.9465369622   * x7
           + 12.9479867399   * x6
           -  4.0507918331   * x5
           -  1.7920253409   * x4
           +  1.4400736722   * x3
           -  0.3545513301   * x2
           +  0.0378074188   * x_norm
           -  0.0014933942;

  // Guard against near-zero denominator
  if (fabs(Q) < 1.0e-30) return D_CB_floor;

  double log10D = P / Q;

  // Clamp log10(D) to physically reasonable range [1e-20, 1e-8] m²/s
  if (log10D < -20.0) log10D = -20.0;
  if (log10D >  -8.0) log10D =  -8.0;

  return pow(10.0, log10D);
}

/* ---------------------------------------------------------------------- */
/*  Pair diffusion coefficient: limiting species rule.
 *    D_pair(i,j) = min( D_particle(i), D_particle(j) )
 * ---------------------------------------------------------------------- */

double FixLithiumDiffusionImplicit::get_diffusion_coefficient(int i, int j)
{
  double Di = compute_D_particle(i);
  double Dj = compute_D_particle(j);
  return std::min(Di, Dj);
}

/* ---------------------------------------------------------------------- */
/*  Compute the chemo-mechanical effective potential for particle i:
 *
 *    Φ_eff,i = U_OCV,i  +  (σ_h,i * Ω_Li,i) / F
 *
 *  σ_h is the hydrostatic stress (Pa, from fix potential/sor)
 *  Ω_Li is the partial molar volume of Li in the host (m³/mol)
 *  F    is Faraday's constant (C/mol)
 *
 *  Compressive stress (σ_h < 0) → Lowers Φ_eff → pushes Li away
 *  Tensile stress    (σ_h > 0) → Raises Φ_eff → attracts Li
 * ---------------------------------------------------------------------- */

double FixLithiumDiffusionImplicit::compute_effective_potential(int i)
{
  double U_ocv = equilibrium_potential[i];
  double sigma_h = hydrostatic_stress[i];
  double Omega = get_Omega_Li(atom->type[i]);

  double effective_potential = U_ocv + (sigma_h * Omega) / F;

  if (effective_potential < 0.0) {
    return 0.0;
  }
  return effective_potential;
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusionImplicit::implicit_diffusion_step(double dt_diff)
{
  int i, j, ii, jj, inum, jnum;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq, r;

  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // ------------------------------------------------------------------
  // STEP 1: Save current concentration as c_old (the c^n term)
  // ------------------------------------------------------------------
  for (i = 0; i < nlocal + atom->nghost; i++) {
    c_old[i] = lithium_concentration[i];
  }

  // ------------------------------------------------------------------
  // STEP 2: Compute cross-type flux (EXPLICIT in chemo-mechanical potential)
  // ------------------------------------------------------------------
  for (i = 0; i < nlocal; i++) {
    cross_flux_arr[i] = 0.0;
    lithium_flux[i] = 0.0;
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    int type_i = type[i];

    if (type_i != AM_type && type_i != CB_type) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    double c_i = lithium_concentration[i];
    double Phi_eff_i = compute_effective_potential(i);

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      int type_j = type[j];
      if (type_j != AM_type && type_j != CB_type && type_j != LM_type) continue;

      // Only cross-type pairs here
      if (type_j == type_i) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      r = sqrt(rsq);

      if (r < (radius[i] + radius[j])) {
        double A = calculate_contact_area(i, j);
        double r_SI = r * 1.0e-6;

        double c_j = (type_j == LM_type) ? c_li_max_LM : lithium_concentration[j];
        double D_ij = get_diffusion_coefficient(i, j);

        double Phi_eff_j = compute_effective_potential(j);

        double M_ij = compute_mobility(c_i, c_j, D_ij);
        double flux = A * M_ij * (Phi_eff_i - Phi_eff_j) / r_SI;

        cross_flux_arr[i] += flux;
      }
    }

    // Subtract electrolyte current contribution
    cross_flux_arr[i] -= current_el[i] / F;
  }

  // ------------------------------------------------------------------
  // STEP 3: Jacobi iteration for implicit same-type diffusion
  // ------------------------------------------------------------------

  int iter;
  double global_max_change;

  for (iter = 0; iter < max_iter; iter++) {
    double local_max_change = 0.0;

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      int type_i = type[i];

      if (type_i != AM_type && type_i != CB_type) continue;

      double V_eff = initial_volume[i] * get_V_exp_max(type_i);
      if (V_eff <= 0.0) continue;

      double diag = 1.0;
      double rhs = c_old[i] + dt_diff * cross_flux_arr[i] / V_eff;

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];

      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        int type_j = type[j];

        // Only same-type neighbors for implicit treatment
        if (type_j != type_i) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        r = sqrt(rsq);

        if (r < (radius[i] + radius[j])) {
          double A = calculate_contact_area(i, j);
          double r_SI = r * 1.0e-6;

          // SOC-dependent pair coefficient (min of the two particles)
          double D_same = get_diffusion_coefficient(i, j);

          double alpha = dt_diff * A * D_same / (r_SI * V_eff);

          diag += alpha;
          rhs += alpha * lithium_concentration[j];
        }
      }

      c_new[i] = rhs / diag;

      double change = fabs(c_new[i] - lithium_concentration[i]);
      if (change > local_max_change) local_max_change = change;
    }

    // Update local concentrations
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if (type[i] == AM_type || type[i] == CB_type) {
        lithium_concentration[i] = c_new[i];
      }
    }

    fix_lithium_concentration->do_forward_comm();

    MPI_Allreduce(&local_max_change, &global_max_change, 1,
                  MPI_DOUBLE, MPI_MAX, world);

    if (global_max_change < tol) {
      iter++;
      break;
    }
  }

  last_iter_count = iter;
  last_residual = global_max_change;

  // ------------------------------------------------------------------
  // STEP 4: Update derived quantities from converged concentration
  // ------------------------------------------------------------------
  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    int type_i = type[i];

    if (type_i == LM_type) {
      lithium_concentration[i] = c_li_max_LM;
      lithium_content[i] = max_lithium_content;
      continue;
    }

    if (type_i == AM_type || type_i == CB_type) {
      double c_max = get_c_li_max(type_i);
      double V_exp_max_cur = get_V_exp_max(type_i);

      lithium_content[i] = lithium_concentration[i] * max_lithium_content / c_max;

      if (lithium_content[i] < initial_lithium_content)
        lithium_content[i] = initial_lithium_content;
      if (lithium_content[i] > target_lithium_content)
        lithium_content[i] = target_lithium_content;

      lithium_concentration[i] = lithium_content[i] * c_max / max_lithium_content;

      li_mols[i] = lithium_content[i] * c_max * initial_volume[i] *
                   V_exp_max_cur / max_lithium_content;

      double V_eff = initial_volume[i] * V_exp_max_cur;
      lithium_flux[i] = (lithium_concentration[i] - c_old[i]) * V_eff / dt_diff;
    }
  }

  // Refresh per-atom D after lithium content changed
  update_per_atom_D();
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusionImplicit::ensure_arrays_sized()
{
  int needed = atom->nmax;
  if (needed > nmax_alloc) {
    memory->grow(c_old, needed, "diff_impl:c_old");
    memory->grow(c_new, needed, "diff_impl:c_new");
    memory->grow(cross_flux_arr, needed, "diff_impl:cross_flux");
    nmax_alloc = needed;
  }
}

/* ---------------------------------------------------------------------- */

void FixLithiumDiffusionImplicit::updatePtrs()
{
  lithium_content = fix_lithium_content->vector_atom;
  lithium_concentration = fix_lithium_concentration->vector_atom;
  equilibrium_potential = fix_equilibrium_potential->vector_atom;
  hydrostatic_stress = fix_hydrostatic_stress->vector_atom;
  current_el = fix_current_el->vector_atom;
  diffusion_coefficient = fix_diffusion_coefficient->vector_atom;
  lithium_flux = fix_lithium_flux->vector_atom;
  li_mols = fix_li_mols->vector_atom;
  initial_volume = fix_initial_volume->vector_atom;
}

/* ---------------------------------------------------------------------- */

double FixLithiumDiffusionImplicit::get_c_li_max(int t)
{
  if (t == AM_type) return c_li_max_AM;
  if (t == CB_type) return c_li_max_CB;
  if (t == LM_type) return c_li_max_LM;
  return 0.0;
}

/* ---------------------------------------------------------------------- */

double FixLithiumDiffusionImplicit::get_V_exp_max(int t)
{
  if (t == AM_type) return V_exp_max_AM;
  if (t == CB_type) return V_exp_max_CB;
  if (t == LM_type) return V_exp_max_LM;
  return 0.0;
}

/* ---------------------------------------------------------------------- */

double FixLithiumDiffusionImplicit::get_Omega_Li(int t)
{
  if (t == AM_type) return Omega_Li_AM;
  if (t == CB_type) return Omega_Li_CB;
  if (t == LM_type) return Omega_Li_LM;
  return 0.0;
}

/* ---------------------------------------------------------------------- */

double FixLithiumDiffusionImplicit::compute_mobility(double c_i, double c_j, double D_ij)
{
  double c_avg = 0.5 * (c_i + c_j);
  return (c_avg * D_ij * F) / (R * T);
}

/* ---------------------------------------------------------------------- */

double FixLithiumDiffusionImplicit::calculate_contact_area(int i, int j)
{
  double *radius = atom->radius;
  double **x = atom->x;
  double L = 1.0e-6;
  double dx = (x[i][0] - x[j][0]) * L;
  double dy = (x[i][1] - x[j][1]) * L;
  double dz = (x[i][2] - x[j][2]) * L;
  double r = sqrt(dx*dx + dy*dy + dz*dz);
  double ri = radius[i] * L;
  double rj = radius[j] * L;
  double rsum = ri + rj;

  if (r >= rsum) return 0.0;

  double delta = rsum - r;
  double reff = ri * rj / rsum;
  return M_PI * delta * reff;
}

/* ---------------------------------------------------------------------- */

double FixLithiumDiffusionImplicit::compute_scalar()
{
  return total_diffusion_time;
}

/* ---------------------------------------------------------------------- */

double FixLithiumDiffusionImplicit::compute_vector(int n)
{
  if (n == 0) return static_cast<double>(last_iter_count);
  if (n == 1) return last_residual;
  if (n == 2) return total_diffusion_time;
  return 0.0;
}