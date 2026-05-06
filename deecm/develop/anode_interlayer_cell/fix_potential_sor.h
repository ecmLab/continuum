/* ----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine
    DCS Computing GmbH, Linz, Austria
    http://www.dcs-computing.com

    Contributing author and copyright for this file:
    Created by: Joseph Vazquez
    Copyright 2024-     DCS Computing GmbH, Linz

    Notes:
    - Solves BOTH electronic (phi_ed) and electrolyte (phi_el) potentials
      via Successive Over-Relaxation (SOR).
    - Boundary conditions at each interface (Interlayer–LM and
      Interlayer–SE) can independently be Dirichlet (fixed potential)
      or Neumann (fixed current density).
    - Current conservation is enforced: total current I = i_app * area.

    -----------------------------------------------------------------------
    Example input-script usage (all keyword/value pairs are optional;
    defaults shown):

      fix pot all potential/sor  &
          omega       1.5                          &
          max_iter    50                           &
          temperature 303.0                        &
          AM_type 1  CB_type 2  LM_type 3  SE_type 4  &
          area        4.0e-10                      &
          i_app       1.0                          &
          sigma_el    0.0  0.0  0.0  0.05          &
          sigma_ed    1.0e5  100.0  1.0e7  0.0     &
          BC_LM_el    potential  0.0                &
          BC_LM_ed    potential  0.0                &
          BC_SE_el    current    1.0                &
          BC_SE_ed    current    1.0

    Keyword reference
    -----------------
    omega          SOR relaxation factor (0 < omega < 2)
    max_iter       SOR iterations per timestep
    temperature    Kelvin
    AM_type        atom type for active-material particles
    CB_type        atom type for carbon-binder particles
    LM_type        atom type for lithium-metal particles
    SE_type        atom type for solid-electrolyte particles
    area           cross-sectional cell area [m^2]
    i_app          applied current density [A/m^2]
    sigma_el       ionic conductivity per type: AM CB LM SE [S/m]
    sigma_ed       electronic conductivity per type: AM CB LM SE [S/m]
    BC_LM_el       electrolyte BC at interlayer–LM interface
    BC_LM_ed       electronic  BC at interlayer–LM interface
    BC_SE_el       electrolyte BC at interlayer–SE interface
    BC_SE_ed       electronic  BC at interlayer–SE interface
                   Each takes "potential <V>" or "current <A/m^2>"
                   (for current BCs the value multiplies i_app)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(potential/sor,FixPotentialSOR)
// clang-format on
#else

#ifndef LMP_FIX_POTENTIAL_SOR_H
#define LMP_FIX_POTENTIAL_SOR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPotentialSOR : public Fix {
 public:
  FixPotentialSOR(class LAMMPS *, int, char **);
  virtual ~FixPotentialSOR();
  virtual void post_create();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void pre_force(int);
  void post_force(int);
  double compute_scalar();
  double compute_vector(int);

  virtual void updatePtrs();

 protected:

  /* ---------- enumerations ------------------------------------------ */
  enum BCMode { BC_POTENTIAL = 0, BC_CURRENT = 1 };

  /* ---------- SOR parameters ---------------------------------------- */
  double omega;
  int    max_iterations;
  int    current_iteration;
  double convergence_el;
  double convergence_ed;

  /* ---------- physical constants ------------------------------------- */
  double R;   // gas constant  [J/(mol·K)]
  double T;   // temperature   [K]
  double FF;  // Faraday const [C/mol]

  /* ---------- particle types ---------------------------------------- */
  int AM_type;   // active material   (default 1)
  int CB_type;   // carbon / binder   (default 2)
  int LM_type;   // lithium metal     (default 3)
  int SE_type;   // solid electrolyte (default 4)

  /* ---------- cell geometry ----------------------------------------- */
  double area;   // cross-sectional area [m^2]
  double i_app;  // applied current density [A/m^2]

  /* ---------- per-type conductivities (indexed 0..3 = AM,CB,LM,SE) -- */
  static const int NTYPES_MAX = 4;
  double sigma_el_type[NTYPES_MAX];  // ionic   [S/m]
  double sigma_ed_type[NTYPES_MAX];  // electronic [S/m]

  /* ---------- boundary-condition specification ----------------------- */
  //   _LM = interlayer – lithium-metal interface
  //   _SE = interlayer – solid-electrolyte interface
  BCMode bc_LM_el_mode;   double bc_LM_el_value;
  BCMode bc_LM_ed_mode;   double bc_LM_ed_value;
  BCMode bc_SE_el_mode;   double bc_SE_el_value;
  BCMode bc_SE_ed_mode;   double bc_SE_ed_value;

  /* ---------- current distribution ---------------------------------- */
  double total_current;
  double global_area_LM;
  double global_area_SE;
  // per-interface, per-potential current densities (computed from i_app)
  double i_dens_LM_el, i_dens_LM_ed;
  double i_dens_SE_el, i_dens_SE_ed;

  /* ---------- per-atom property pointers ----------------------------- */
  double *phi_el,     *phi_el_old;
  double *phi_ed,     *phi_ed_old;
  double *current_el, *current_ed;
  double *hydrostatic_stress;
  double *contact_area_LM, *contact_area_SE;

  /* ---------- fix property/atom handles ------------------------------ */
  class FixPropertyAtom *fix_phi_el,     *fix_phi_el_old;
  class FixPropertyAtom *fix_phi_ed,     *fix_phi_ed_old;
  class FixPropertyAtom *fix_current_el, *fix_current_ed;
  class FixPropertyAtom *fix_hydrostatic_stress;
  class FixPropertyAtom *fix_init_flag;
  class FixPropertyAtom *fix_contact_area_LM;
  class FixPropertyAtom *fix_contact_area_SE;

  /* ---------- neighbour list ---------------------------------------- */
  class NeighList *list;
  bool first_run;

  /* ---------- internal helpers -------------------------------------- */

  /// Map an atom type (1-based) to local type-index (0-based) or -1.
  inline int type_index(int t) const {
    if (t == AM_type) return 0;
    if (t == CB_type) return 1;
    if (t == LM_type) return 2;
    if (t == SE_type) return 3;
    return -1;
  }
  inline bool is_interlayer(int t) const { return (t == AM_type || t == CB_type); }
  inline bool is_LM(int t) const { return (t == LM_type); }
  inline bool is_SE(int t) const { return (t == SE_type); }
  inline bool is_LM_interface(int ti, int tj) const {
    return (is_interlayer(ti) && is_LM(tj)) || (is_LM(ti) && is_interlayer(tj));
  }
  inline bool is_SE_interface(int ti, int tj) const {
    return (is_interlayer(ti) && is_SE(tj)) || (is_SE(ti) && is_interlayer(tj));
  }

  /// Harmonic mean of two conductivities (returns 0 if either is 0).
  inline double harmonic_mean(double a, double b) const {
    if (a < 1e-30 || b < 1e-30) return 0.0;
    return 2.0 * a * b / (a + b);
  }

  void calculate_interface_areas();
  void solve_potential_iteration(double *phi, double *phi_old,
                                 const double *sigma_type,
                                 BCMode bc_LM_mode, double bc_LM_val,
                                 BCMode bc_SE_mode, double bc_SE_val,
                                 double i_dens_LM, double i_dens_SE,
                                 double *current_arr,
                                 FixPropertyAtom *fix_phi);
  void apply_boundary_conditions();
  void calculate_hydrostatic_stress();
  double calculate_contact_area(int, int);
  double check_convergence(double *phi, double *phi_old);
};

}  // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix potential/sor command
  Check syntax against the keyword reference above.

E: Fix potential/sor requires newton pair off
  This fix requires newton pair off.

E: Could not find required property/atom fixes
  Internal error – missing property/atom fixes.

E: Invalid particle types for potential/sor
  AM_type, CB_type, LM_type, SE_type must be valid atom types.

E: Current BCs on both sides without a Dirichlet reference
  At least one interface must use a potential (Dirichlet) BC for each
  field to anchor the solution.

E: FixPotentialSOR: Could not find compute with ID 'st'
  Add "compute st all stress/atom" before this fix.
*/