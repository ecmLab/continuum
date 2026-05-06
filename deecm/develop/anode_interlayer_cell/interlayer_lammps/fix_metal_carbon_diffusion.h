/* ----------------------------------------------------------------------
    LIGGGHTS® - DEM simulation engine
    http://www.dcs-computing.com, office@dcs-computing.com
    LIGGGHTS® is part of CFDEM®project: http://www.liggghts.com | http://www.cfdem.com

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    Created: Joseph Vazquez Mercado, RIT 2025
    Copyright 2024-     DCS Computing GmbH, Linz
    Notes: Metal-Carbon system lithium diffusion header
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(metal_carbon_diffusion,FixMetalCarbonDiffusion)

#else

#ifndef LMP_FIX_METAL_CARBON_DIFFUSION_H
#define LMP_FIX_METAL_CARBON_DIFFUSION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPropertyAtom;
class NeighList;

class FixMetalCarbonDiffusion : public Fix {
 public:
  FixMetalCarbonDiffusion(class LAMMPS *, int, char **);
  virtual ~FixMetalCarbonDiffusion();
  void post_create();
  int setmask();
  void init();
  void setup(int);
  void init_list(int, class NeighList *);
  void pre_force(int);
  void post_force(int);
  double compute_scalar();
  void updatePtrs();

 protected:
  // Diffusion parameters
  double c_li_max;           // Maximum Li concentration in carbon (mol/m³)
  double c_li_source;        // Li concentration at source region (mol/m³)
  double D_Li_C;             // Li diffusion coefficient in carbon (m²/s)
  double D_Li_metal;         // Li diffusion coefficient in metal (m²/s)

  // Source region parameters
  enum SourceType {
    SOURCE_NONE = 0,
    SOURCE_ZLO,
    SOURCE_ZHI,
    SOURCE_XLO,
    SOURCE_XHI,
    SOURCE_YLO,
    SOURCE_YHI
  };
  int source_type;
  double source_coord;       // Coordinate of source plane (μm)
  double source_thickness;   // Thickness of source region (μm)

  // Property pointers
  double *lithium_content;      // Li/Metal molar ratio
  double *lithium_concentration; // mol/m³
  double *diffusion_coefficient; // m²/s
  double *lithium_flux;         // mol/s
  double *li_mols;              // mols of Li in particle
  double *metal_mols;           // mols of metal in particle
  
  // Fix pointers
  class FixPropertyAtom *fix_lithium_content;
  class FixPropertyAtom *fix_lithium_concentration;
  class FixPropertyAtom *fix_diffusion_coefficient;
  class FixPropertyAtom *fix_lithium_flux;
  class FixPropertyAtom *fix_li_mols;
  class FixPropertyAtom *fix_metal_mols;
  
  // Particle types
  int metal_type;   // Atom type for metal particles (default: 1)
  int carbon_type;  // Atom type for carbon particles (default: 2)
  
  // Neighbor list
  class NeighList *list;
  
  // Methods
  void update_lithium_content();
  double calculate_contact_area(int, int);
  bool is_in_source_region(double *pos);
  double distance_to_source(double *pos);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix metal_carbon_diffusion command

Self-explanatory. Check the input script syntax and compare to the
documentation for the command.

E: Fix metal_carbon_diffusion requires lithiumContent property

The lithiumContent property/atom must be available.

E: Fix metal_carbon_diffusion requires lithiumConcentration property

The lithiumConcentration property/atom must be available.

E: Could not find required property/atom fixes

Required property fixes (diffusionCoefficient, lithiumFlux, lithiumMols, metalMols)
were not found.

E: Unknown source type for fix metal_carbon_diffusion

Valid source types are: zlo, zhi, xlo, xhi, ylo, yhi

*/