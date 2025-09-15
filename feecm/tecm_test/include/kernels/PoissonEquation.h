#pragma once

#include "ADKernel.h"

/**
 * PoissonEquation implements the space charge Poisson equation for TECM-no-EN
 * Based on the derivation in TECM_noEN_transference.pdf
 * 
 * Poisson equation:
 * -∇ · (ε ∇Ψ) = -ρₑ
 * 
 * Where:
 * - ε is the permittivity 
 * - Ψ is the Poisson potential
 * - ρₑ is the charge density = F(z₊c₊ + z₋c₋) + ρfixed
 */
class PoissonEquation : public ADKernel
{
public:
  static InputParameters validParams();

  PoissonEquation(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// Permittivity material property
  const ADMaterialProperty<Real> & _permittivity;

  /// Charge density variable (ρₑ)
  const ADVariableValue & _charge_density;
};