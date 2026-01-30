#pragma once
#include "ADKernel.h"

/// -∫ w * rho_e dOmega, with rho_e = F (z_plus c_plus + z_minus c_minus)
class ADPoissonChargeRHS : public ADKernel
{
public:
  static InputParameters validParams();
  ADPoissonChargeRHS(const InputParameters & parameters);

protected:
  ADReal computeQpResidual() override;

  const ADVariableValue & _c_plus;
  const ADVariableValue & _c_minus;

  const Real _F;      // Faraday constant
  const Real _z_plus; // +1
  const Real _z_minus;// -1
};
