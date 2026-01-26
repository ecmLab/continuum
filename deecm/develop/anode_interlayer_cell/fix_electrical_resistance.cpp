/* ----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine
    
    Implements the electro-mechanical contact model from:
    Zhang et al. (2024) "An electro-mechanical contact model for 
    particulate systems", Powder Technology 440, 119759
    
    Total resistance R_ij = R_i + R_c,ij + R_j
    where R_i, R_j are bulk resistances and R_c,ij is contact resistance
    
    Full MPI support with ghost atom communication.
    
    UNITS:
    - LIGGGHTS micro system: lengths in μm
    - All internal calculations: SI units (m, Ω, S)
    - Input conductivity: S/m
    - Output resistance: Ω
------------------------------------------------------------------------- */

#include "fix_electrical_resistance.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "error.h"
#include "force.h"
#include "pair_gran.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "fix_property_atom.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include <cmath>
#include <cstring>
#include <algorithm>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixElectricalResistance::FixElectricalResistance(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  T_(303.0),
  AM_type_(1),
  CB_type_(2),
  SE_type_(3),
  LM_type_(4),
  wall_type_(5),
  ntypes_(5),
  wall_resistivity_(1.0e-7),      // Default: ~10^7 S/m (metallic)
  min_overlap_(1.0e-12),          // Default: 1e-6 μm = 1e-12 m
  total_conductance_(nullptr),
  total_resistance_(nullptr),
  contact_count_(nullptr),
  avg_contact_resistance_(nullptr),
  fix_total_conductance_(nullptr),
  fix_total_resistance_(nullptr),
  fix_contact_count_(nullptr),
  fix_avg_contact_resistance_(nullptr),
  pair_gran_(nullptr),
  list_(nullptr),
  newton_pair_(0),
  nmax_comm_(0),
  total_system_resistance_(0.0),
  avg_pair_resistance_(0.0),
  total_contact_pairs_(0)
{
  if (narg < 3)
    error->all(FLERR, "Illegal fix electrical_resistance command - too few arguments");

  // Initialize resistivity vector with default values [Ω·m]
  // Default: σ = 10^7 S/m → ρ = 10^-7 Ω·m (metallic conductor)
  resistivity_.resize(ntypes_ + 1, 1.0e-7);

  // Parse arguments
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "temperature") == 0) {
      if (iarg + 2 > narg) 
        error->all(FLERR, "Illegal fix electrical_resistance command - missing temperature value");
      T_ = force->numeric(FLERR, arg[iarg + 1]);
      iarg += 2;
    } 
    else if (strcmp(arg[iarg], "AM_type") == 0) {
      if (iarg + 2 > narg) 
        error->all(FLERR, "Illegal fix electrical_resistance command - missing AM_type value");
      AM_type_ = force->inumeric(FLERR, arg[iarg + 1]);
      iarg += 2;
    } 
    else if (strcmp(arg[iarg], "CB_type") == 0) {
      if (iarg + 2 > narg) 
        error->all(FLERR, "Illegal fix electrical_resistance command - missing CB_type value");
      CB_type_ = force->inumeric(FLERR, arg[iarg + 1]);
      iarg += 2;
    } 
    else if (strcmp(arg[iarg], "SE_type") == 0) {
      if (iarg + 2 > narg) 
        error->all(FLERR, "Illegal fix electrical_resistance command - missing SE_type value");
      SE_type_ = force->inumeric(FLERR, arg[iarg + 1]);
      iarg += 2;
    } 
    else if (strcmp(arg[iarg], "LM_type") == 0) {
      if (iarg + 2 > narg) 
        error->all(FLERR, "Illegal fix electrical_resistance command - missing LM_type value");
      LM_type_ = force->inumeric(FLERR, arg[iarg + 1]);
      iarg += 2;
    }
    else if (strcmp(arg[iarg], "wall_type") == 0) {
      if (iarg + 2 > narg) 
        error->all(FLERR, "Illegal fix electrical_resistance command - missing wall_type value");
      wall_type_ = force->inumeric(FLERR, arg[iarg + 1]);
      iarg += 2;
    }
    else if (strcmp(arg[iarg], "ntypes") == 0) {
      if (iarg + 2 > narg) 
        error->all(FLERR, "Illegal fix electrical_resistance command - missing ntypes value");
      ntypes_ = force->inumeric(FLERR, arg[iarg + 1]);
      resistivity_.resize(ntypes_ + 1, 1.0e-7);
      iarg += 2;
    }
    else if (strcmp(arg[iarg], "conductivity") == 0) {
      // Expects: conductivity type1 sigma1 type2 sigma2 ...
      // sigma in S/m, converted to resistivity in Ω·m: ρ = 1/σ
      iarg++;
      while (iarg + 1 < narg) {
        char first_char = arg[iarg][0];
        if (!isdigit(first_char) && first_char != '-' && first_char != '+')
          break;
        
        int type = force->inumeric(FLERR, arg[iarg]);
        double sigma = force->numeric(FLERR, arg[iarg + 1]);
        
        if (type < 1 || type > ntypes_)
          error->all(FLERR, "Invalid type in conductivity specification");
        if (sigma <= 0.0)
          error->all(FLERR, "Conductivity must be positive");
        
        // Convert conductivity to resistivity: ρ[Ω·m] = 1 / σ[S/m]
        resistivity_[type] = 1.0 / sigma;
        iarg += 2;
      }
    }
    else if (strcmp(arg[iarg], "wall_conductivity") == 0) {
      if (iarg + 2 > narg) 
        error->all(FLERR, "Illegal fix electrical_resistance command - missing wall_conductivity value");
      double sigma_wall = force->numeric(FLERR, arg[iarg + 1]);
      if (sigma_wall <= 0.0)
        error->all(FLERR, "Wall conductivity must be positive");
      // Convert: ρ[Ω·m] = 1 / σ[S/m]
      wall_resistivity_ = 1.0 / sigma_wall;
      iarg += 2;
    }
    else if (strcmp(arg[iarg], "min_overlap") == 0) {
      // Input is in μm (LIGGGHTS micro units), convert to m
      if (iarg + 2 > narg) 
        error->all(FLERR, "Illegal fix electrical_resistance command - missing min_overlap value");
      double min_overlap_um = force->numeric(FLERR, arg[iarg + 1]);
      min_overlap_ = min_overlap_um * LENGTH_CONV;  // Convert μm to m
      iarg += 2;
    }
    else {
      char errstr[256];
      sprintf(errstr, "Unknown keyword '%s' in fix electrical_resistance command", arg[iarg]);
      error->all(FLERR, errstr);
    }
  }

  // Setup compute flags
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 4;
  global_freq = 1;
  extscalar = 0;
  extvector = 0;

  // Enable communication
  comm_forward = 4;
  comm_reverse = 2;
}

/* ---------------------------------------------------------------------- */

FixElectricalResistance::~FixElectricalResistance()
{
}

/* ---------------------------------------------------------------------- */

void FixElectricalResistance::post_create()
{
  // Create per-atom property for total conductance [S]
  fix_total_conductance_ = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("totalConductance", "property/atom", "scalar", 0, 0, style, false));
  if (!fix_total_conductance_) {
    const char* fixarg[9];
    fixarg[0] = "totalConductance";
    fixarg[1] = "all";
    fixarg[2] = "property/atom";
    fixarg[3] = "totalConductance";
    fixarg[4] = "scalar";
    fixarg[5] = "no";
    fixarg[6] = "yes";
    fixarg[7] = "yes";
    fixarg[8] = "0.0";
    fix_total_conductance_ = modify->add_fix_property_atom(9, const_cast<char**>(fixarg), style);
  }

  // Create per-atom property for total resistance [Ω]
  fix_total_resistance_ = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("totalResistance", "property/atom", "scalar", 0, 0, style, false));
  if (!fix_total_resistance_) {
    const char* fixarg[9];
    fixarg[0] = "totalResistance";
    fixarg[1] = "all";
    fixarg[2] = "property/atom";
    fixarg[3] = "totalResistance";
    fixarg[4] = "scalar";
    fixarg[5] = "no";
    fixarg[6] = "yes";
    fixarg[7] = "no";
    fixarg[8] = "1.0e20";
    fix_total_resistance_ = modify->add_fix_property_atom(9, const_cast<char**>(fixarg), style);
  }

  // Create per-atom property for contact count
  fix_contact_count_ = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("resistanceContactCount", "property/atom", "scalar", 0, 0, style, false));
  if (!fix_contact_count_) {
    const char* fixarg[9];
    fixarg[0] = "resistanceContactCount";
    fixarg[1] = "all";
    fixarg[2] = "property/atom";
    fixarg[3] = "resistanceContactCount";
    fixarg[4] = "scalar";
    fixarg[5] = "no";
    fixarg[6] = "yes";
    fixarg[7] = "yes";
    fixarg[8] = "0.0";
    fix_contact_count_ = modify->add_fix_property_atom(9, const_cast<char**>(fixarg), style);
  }

  // Create per-atom property for average contact resistance [Ω]
  fix_avg_contact_resistance_ = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("avgContactResistance", "property/atom", "scalar", 0, 0, style, false));
  if (!fix_avg_contact_resistance_) {
    const char* fixarg[9];
    fixarg[0] = "avgContactResistance";
    fixarg[1] = "all";
    fixarg[2] = "property/atom";
    fixarg[3] = "avgContactResistance";
    fixarg[4] = "scalar";
    fixarg[5] = "no";
    fixarg[6] = "yes";
    fixarg[7] = "no";
    fixarg[8] = "0.0";
    fix_avg_contact_resistance_ = modify->add_fix_property_atom(9, const_cast<char**>(fixarg), style);
  }
}

/* ---------------------------------------------------------------------- */

int FixElectricalResistance::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixElectricalResistance::init()
{
  // Get newton pair setting
  newton_pair_ = force->newton_pair;

  // Get pair style (verify granular simulation)
  pair_gran_ = static_cast<PairGran*>(force->pair_match("gran", 0));
  if (!pair_gran_)
    error->all(FLERR, "Fix electrical_resistance requires a granular pair style");

  // Find per-atom property fixes
  fix_total_conductance_ = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("totalConductance", "property/atom", "scalar", 0, 0, style));
  fix_total_resistance_ = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("totalResistance", "property/atom", "scalar", 0, 0, style));
  fix_contact_count_ = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("resistanceContactCount", "property/atom", "scalar", 0, 0, style));
  fix_avg_contact_resistance_ = static_cast<FixPropertyAtom*>(
    modify->find_fix_property("avgContactResistance", "property/atom", "scalar", 0, 0, style));

  // Request neighbor list - full list with ghost atoms
  irequest_ = neighbor->request(this);
  neighbor->requests[irequest_]->pair = 0;
  neighbor->requests[irequest_]->fix = 1;
  neighbor->requests[irequest_]->half = 0;
  neighbor->requests[irequest_]->full = 1;
  neighbor->requests[irequest_]->occasional = 0;
  neighbor->requests[irequest_]->ghost = 1;

  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixElectricalResistance::init_list(int id, NeighList *ptr)
{
  list_ = ptr;
}

/* ---------------------------------------------------------------------- */

void FixElectricalResistance::setup(int vflag)
{
  updatePtrs();
  
  if (fix_total_conductance_) fix_total_conductance_->do_forward_comm();
  if (fix_total_resistance_) fix_total_resistance_->do_forward_comm();
  if (fix_contact_count_) fix_contact_count_->do_forward_comm();
  if (fix_avg_contact_resistance_) fix_avg_contact_resistance_->do_forward_comm();
  
  calculate_resistances();
  
  if (fix_total_conductance_) fix_total_conductance_->do_forward_comm();
  if (fix_total_resistance_) fix_total_resistance_->do_forward_comm();
  if (fix_contact_count_) fix_contact_count_->do_forward_comm();
  if (fix_avg_contact_resistance_) fix_avg_contact_resistance_->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

void FixElectricalResistance::pre_force(int vflag)
{
  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixElectricalResistance::post_force(int vflag)
{
  calculate_resistances();
  
  if (fix_total_conductance_) fix_total_conductance_->do_forward_comm();
  if (fix_total_resistance_) fix_total_resistance_->do_forward_comm();
  if (fix_contact_count_) fix_contact_count_->do_forward_comm();
  if (fix_avg_contact_resistance_) fix_avg_contact_resistance_->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

void FixElectricalResistance::updatePtrs()
{
  if (fix_total_conductance_)
    total_conductance_ = fix_total_conductance_->vector_atom;
  if (fix_total_resistance_)
    total_resistance_ = fix_total_resistance_->vector_atom;
  if (fix_contact_count_)
    contact_count_ = fix_contact_count_->vector_atom;
  if (fix_avg_contact_resistance_)
    avg_contact_resistance_ = fix_avg_contact_resistance_->vector_atom;
}

/* ---------------------------------------------------------------------- */

double FixElectricalResistance::get_resistivity(int type)
{
  // Returns resistivity in SI units [Ω·m]
  if (type < 1 || type > ntypes_) return 1.0e20;
  return resistivity_[type];
}

/* ---------------------------------------------------------------------- */

double FixElectricalResistance::get_effective_radius(double r_i, double r_j)
{
  // r_eff = r_i * r_j / (r_i + r_j)
  // Input and output in SI units [m]
  return (r_i * r_j) / (r_i + r_j);
}

/* ---------------------------------------------------------------------- */

double FixElectricalResistance::compute_overlap(int i, int j, double &L_ij)
{
  // Calculate overlap and center-to-center distance
  // LIGGGHTS provides positions and radii in μm, convert to m
  
  double **x = atom->x;
  double *radius = atom->radius;
  
  // Convert positions from μm to m
  double dx = (x[j][0] - x[i][0]) * LENGTH_CONV;
  double dy = (x[j][1] - x[i][1]) * LENGTH_CONV;
  double dz = (x[j][2] - x[i][2]) * LENGTH_CONV;
  
  // Center-to-center distance [m]
  L_ij = std::sqrt(dx * dx + dy * dy + dz * dz);
  
  // Convert radii from μm to m
  double r_i = radius[i] * LENGTH_CONV;
  double r_j = radius[j] * LENGTH_CONV;
  double radsum = r_i + r_j;
  
  // Overlap (δ) in [m] - positive when particles are in contact
  return radsum - L_ij;
}

/* ---------------------------------------------------------------------- */

double FixElectricalResistance::compute_contact_radius(double overlap, double r_i, double r_j)
{
  // Direct geometric calculation of contact radius
  // All inputs and output in SI units [m]
  //
  // For two overlapping spheres:
  // r_c² = δ × r_eff - δ²/4
  // where δ = overlap, r_eff = r_i × r_j / (r_i + r_j)
  
  if (overlap <= 0.0) return 0.0;
  
  double r_eff = get_effective_radius(r_i, r_j);
  
  // Contact radius squared
  double r_c_sq = overlap * r_eff - 0.25 * overlap * overlap;
  
  if (r_c_sq <= 0.0) return 0.0;
  
  return std::sqrt(r_c_sq);  // [m]
}

/* ---------------------------------------------------------------------- */

double FixElectricalResistance::compute_contact_area(double overlap, double r_i, double r_j)
{
  // Direct calculation of contact area
  // A_c = π × δ × r_eff
  // All inputs in [m], output in [m²]
  
  if (overlap <= 0.0) return 0.0;
  
  double r_eff = get_effective_radius(r_i, r_j);
  return PI * overlap * r_eff;  // [m²]
}

/* ---------------------------------------------------------------------- */

double FixElectricalResistance::compute_contact_resistance(double rho_i, double rho_j, double r_c)
{
  // Holm contact resistance formula (Equation 12 from paper)
  // R_c,ij = (ρ_i + ρ_j) / (4 × r_c)
  //
  // Units: [Ω·m] / [m] = [Ω]
  
  if (r_c <= 0.0) return 1.0e20;
  
  return (rho_i + rho_j) / (4.0 * r_c);  // [Ω]
}

/* ---------------------------------------------------------------------- */

double FixElectricalResistance::compute_bulk_resistance(double rho, double r_i, 
                                                        double r_j, double L_ij)
{
  // Bulk resistance from particle center to contact plane (Equation 17a)
  // R_i = (ρ / (2πr_i)) × ln((r_j² - r_i² - L² + 4r_i×L) / (r_i² - r_j² + L²))
  //
  // All inputs in SI: rho [Ω·m], r_i, r_j, L_ij [m]
  // Output: [Ω]
  
  if (r_i <= 0.0 || L_ij <= 0.0) return 1.0e20;
  
  double r_i_sq = r_i * r_i;
  double r_j_sq = r_j * r_j;
  double L_ij_sq = L_ij * L_ij;
  
  double numerator = r_j_sq - r_i_sq - L_ij_sq + 4.0 * r_i * L_ij;
  double denominator = r_i_sq - r_j_sq + L_ij_sq;
  
  // Check for valid ratio
  if (numerator <= 0.0 || denominator <= 0.0) {
    // Fallback: cylindrical approximation R = ρ × L / (π × r²)
    // Units: [Ω·m] × [m] / [m²] = [Ω]
    return rho * L_ij / (PI * r_i_sq);
  }
  
  double ratio = numerator / denominator;
  if (ratio <= 0.0) return 1.0e20;
  
  // Units: [Ω·m] / [m] = [Ω]
  return (rho / (2.0 * PI * r_i)) * std::log(ratio);  // [Ω]
}

/* ---------------------------------------------------------------------- */

double FixElectricalResistance::compute_total_resistance(int i, int j)
{
  // Total resistance: R_ij = R_i + R_c,ij + R_j (Equation 18)
  // All calculations in SI units, output in [Ω]
  
  double *radius = atom->radius;
  int *type = atom->type;
  
  // Convert radii from μm to m
  double r_i = radius[i] * LENGTH_CONV;
  double r_j = radius[j] * LENGTH_CONV;
  int type_i = type[i];
  int type_j = type[j];
  
  // Get overlap and center-to-center distance [m]
  double L_ij;
  double overlap = compute_overlap(i, j, L_ij);
  
  // Check for valid contact (min_overlap_ is already in [m])
  if (overlap < min_overlap_) return 1.0e20;
  
  // Get material resistivities [Ω·m]
  double rho_i = get_resistivity(type_i);
  double rho_j = get_resistivity(type_j);
  
  // Calculate contact radius from geometry [m]
  double r_c = compute_contact_radius(overlap, r_i, r_j);
  
  // Calculate resistances [Ω]
  // Bulk resistances (Equation 17a, 17b)
  double R_i = compute_bulk_resistance(rho_i, r_i, r_j, L_ij);
  double R_j = compute_bulk_resistance(rho_j, r_j, r_i, L_ij);
  
  // Contact resistance using Holm formula (Equation 12)
  double R_c = compute_contact_resistance(rho_i, rho_j, r_c);
  
  return R_i + R_c + R_j;  // [Ω]
}

/* ---------------------------------------------------------------------- */

void FixElectricalResistance::calculate_resistances()
{
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int *mask = atom->mask;
  double *radius = atom->radius;
  tagint *tag = atom->tag;

  // Reset per-atom values for local AND ghost atoms
  for (int i = 0; i < nall; i++) {
    total_conductance_[i] = 0.0;
    total_resistance_[i] = 1.0e20;
    contact_count_[i] = 0.0;
    avg_contact_resistance_[i] = 0.0;
  }

  // Local accumulators for resistance sum
  std::vector<double> resistance_sum(nall, 0.0);

  // Reset global counters
  int local_contact_pairs = 0;
  double local_sum_pair_resistance = 0.0;
  double local_sum_conductance = 0.0;

  // Get neighbor list
  if (!list_) return;
  
  int inum = list_->inum;
  int *ilist = list_->ilist;
  int *numneigh = list_->numneigh;
  int **firstneigh = list_->firstneigh;

  // Loop over neighbor list
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;

      if (!(mask[j] & groupbit)) continue;

      // Check for contact using direct geometry
      double L_ij;
      double overlap = compute_overlap(i, j, L_ij);
      
      // min_overlap_ is already in [m]
      if (overlap < min_overlap_) continue;

      // Calculate total resistance [Ω]
      double R_ij = compute_total_resistance(i, j);

      if (R_ij > 0.0 && R_ij < 1.0e19) { // if (R_ij > 0.0 && R_ij < 1.0e19) {
        double C_ij = 1.0 / R_ij;  // Conductance [S]
        
        // Accumulate for particle i
        total_conductance_[i] += C_ij;
        resistance_sum[i] += R_ij;
        contact_count_[i] += 1.0;

        // Unique pair counting using tags
        if (i < nlocal) {
          if (j < nlocal) {
            if (tag[i] < tag[j]) {
              local_contact_pairs++;
              local_sum_pair_resistance += R_ij;
              local_sum_conductance += C_ij;
            }
          } else {
            if (tag[i] < tag[j]) {
              local_contact_pairs++;
              local_sum_pair_resistance += R_ij;
              local_sum_conductance += C_ij;
            }
          }
        }
      }
    }
  }

  // Calculate derived per-atom quantities for local atoms
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    
    if (contact_count_[i] > 0) {
      total_resistance_[i] = 1.0 / total_conductance_[i];  // [Ω]
      avg_contact_resistance_[i] = resistance_sum[i] / contact_count_[i];  // [Ω]
    }
  }

  // MPI reduction for global values
  MPI_Allreduce(&local_contact_pairs, &total_contact_pairs_, 1, MPI_INT, MPI_SUM, world);
  
  double global_sum_resistance, global_sum_conductance;
  MPI_Allreduce(&local_sum_pair_resistance, &global_sum_resistance, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&local_sum_conductance, &global_sum_conductance, 1, MPI_DOUBLE, MPI_SUM, world);

  if (total_contact_pairs_ > 0) {
    avg_pair_resistance_ = global_sum_resistance / total_contact_pairs_;  // [Ω]
    total_system_resistance_ = 1.0 / global_sum_conductance;  // [Ω]
  } else {
    avg_pair_resistance_ = 1.0e20;
    total_system_resistance_ = 1.0e20;
  }
}

/* ---------------------------------------------------------------------- */

int FixElectricalResistance::pack_forward_comm(int n, int *list, double *buf,
                                                int pbc_flag, int *pbc)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    buf[m++] = total_conductance_[j];
    buf[m++] = total_resistance_[j];
    buf[m++] = contact_count_[j];
    buf[m++] = avg_contact_resistance_[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixElectricalResistance::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {
    total_conductance_[i] = buf[m++];
    total_resistance_[i] = buf[m++];
    contact_count_[i] = buf[m++];
    avg_contact_resistance_[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int FixElectricalResistance::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {
    buf[m++] = total_conductance_[i];
    buf[m++] = contact_count_[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixElectricalResistance::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    total_conductance_[j] += buf[m++];
    contact_count_[j] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

double FixElectricalResistance::compute_scalar()
{
  // Returns total system resistance [Ω]
  return total_system_resistance_;
}

/* ---------------------------------------------------------------------- */

double FixElectricalResistance::compute_vector(int n)
{
  // All outputs in SI units [Ω] for resistance, [S] for conductance
  // 0: total system resistance (parallel combination) [Ω]
  // 1: average pair resistance [Ω]
  // 2: total contact pairs [-]
  // 3: global average contact resistance per particle [Ω]
  
  if (n == 0) return total_system_resistance_;
  if (n == 1) return avg_pair_resistance_;
  if (n == 2) return static_cast<double>(total_contact_pairs_);
  
  if (n == 3) {
    int nlocal = atom->nlocal;
    int *mask = atom->mask;
    double sum = 0.0;
    int count = 0;
    
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      if (contact_count_[i] > 0) {
        sum += avg_contact_resistance_[i];
        count++;
      }
    }
    
    double global_sum;
    int global_count;
    MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(&count, &global_count, 1, MPI_INT, MPI_SUM, world);
    
    if (global_count > 0) return global_sum / global_count;
    return 0.0;
  }
  
  return 0.0;
}