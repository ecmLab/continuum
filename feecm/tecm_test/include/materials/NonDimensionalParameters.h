// Copyright 2025, ToBeDecided, All Rights Reserved
// License: L-GPL 3.0
#pragma once

#include "Material.h"

/**
 * Material class that provides all characteristic scales and dimensionless parameters
 * for the TECM non-dimensionalization framework.
 * 
 * Based on Li-ion battery characteristic scales:
 * - Length: L₀ = Fc₀D₀/j₀ (diffusion length based on electrochemical kinetics)
 * - Time: t₀ = L₀²/D₀ = F²c₀²D₀/j₀²
 * - Concentration: c₀ = Li-ion concentration in Li metal
 * - Potential: φ₀ = RT₀/F (thermal voltage)
 * - Chemical potential: μ₀ = RT₀
 * - Current density: j₀ = exchange current density at Li metal/SE interface
 * - Stress: σ₀ = RT₀/Ω₀ (where Ω₀ is Li molar volume in Li metal)
 */
class NonDimensionalParameters : public Material
{
public:
  static InputParameters validParams();

  NonDimensionalParameters(const InputParameters & parameters);

  void computeQpProperties() override;

protected:
  ///@{ Characteristic scales (dimensional reference values)
  /// Reference concentration [mol/m³]
  const Real _c0;
  /// Reference temperature [K]  
  const Real _T0;
  /// Reference diffusivity [m²/s]
  const Real _D0;
  /// Exchange current density [A/m²]
  const Real _j0;
  /// Li molar volume in Li metal [m³/mol]
  const Real _Omega0;
  /// Faraday constant [C/mol]
  const Real _F;
  /// Gas constant [J/mol/K]
  const Real _R;
  ///@}

  ///@{ Computed characteristic scales
  /// Characteristic length [m]
  MaterialProperty<Real> & _L0;
  /// Characteristic time [s]
  MaterialProperty<Real> & _t0;
  /// Characteristic potential [V]
  MaterialProperty<Real> & _phi0;
  /// Characteristic chemical potential [J/mol]
  MaterialProperty<Real> & _mu0;
  /// Characteristic stress [Pa]
  MaterialProperty<Real> & _sigma0;
  ///@}

  ///@{ Dimensionless parameters
  /// Electrostatic parameter κ = F²c₀L₀²/(εRT₀)
  MaterialProperty<Real> & _kappa;
  /// Mechanical coupling parameter M = Ω₀c₀
  MaterialProperty<Real> & _M_mech;
  ///@}

  ///@{ Optional material properties for local calculations
  /// Local permittivity for electrostatic parameter
  const MaterialProperty<Real> * _epsilon;
  ///@}
};