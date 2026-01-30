// Copyright 2025, CEWLAB, All Rights Reserved
// License: L-GPL 3.0
#pragma once

#include "ADMaterial.h"

/**
 * Material that computes dimensionless diffusivity D̃ = D/D₀
 * where D₀ is the characteristic diffusivity from NonDimensionalParameters
 */
class NonDimensionalDiffusivity : public ADMaterial
{
public:
  static InputParameters validParams();

  NonDimensionalDiffusivity(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Input dimensional diffusivity
  const Real _D_dimensional;

  /// Characteristic diffusivity D₀
  const MaterialProperty<Real> & _D0;

  /// Output dimensionless diffusivity D̃ = D/D₀
  ADMaterialProperty<Real> & _diffusivity_nondim;

  /// Optional: output the dimensional diffusivity for comparison
  ADMaterialProperty<Real> & _diffusivity_dim;
};