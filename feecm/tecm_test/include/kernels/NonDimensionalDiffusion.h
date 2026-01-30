// Copyright 2025, CEWLAB, All Rights Reserved
// License: L-GPL 3.0
#pragma once

#include "ADKernel.h"

/**
 * Non-dimensional diffusion kernel for the TECM framework.
 * 
 * Implements: -∇̃ · (D̃ ∇̃ũ) where:
 * - ũ is the dimensionless variable (concentration, potential, etc.)
 * - D̃ = D/D₀ is the dimensionless diffusivity
 * - ∇̃ = L₀∇ is the dimensionless gradient operator
 * - All quantities are automatically scaled using characteristic scales from NonDimensionalParameters
 *
 * The weak form becomes: (D̃ ∇̃ũ, ∇̃φ) where φ is the test function
 */
class NonDimensionalDiffusion : public ADKernel
{
public:
  static InputParameters validParams();

  NonDimensionalDiffusion(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// Dimensionless diffusivity D̃ = D/D₀
  const ADMaterialProperty<Real> & _diffusivity_nondim;

  /// Characteristic length scale L₀ (for reference, though scaling handled in weak form)
  const MaterialProperty<Real> & _L0;

  /// Whether to use automatic scaling (default true)
  const bool _use_nondimensional_scaling;

  /// Optional dimensional diffusivity for comparison mode
  const ADMaterialProperty<Real> * _diffusivity_dim;
};