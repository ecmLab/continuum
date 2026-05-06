// SternRobinBC: Robin boundary condition for electrostatic potential
// Standard MOOSE BC implementation replaces the natural flux term
//   -\int_Gamma (epsilon * \partial_n u) w dA
// using the Robin relation epsilon * \partial_n u = C_S * (phi_M - u),
// yielding the boundary residual integrand
//   (C_S * (u - phi_M)) * w.

#pragma once

#include "ADIntegratedBC.h"

class SternRobinBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();
  SternRobinBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

private:
  const Real _C_S;   // Computed Stern capacitance per area [F/m^2]
  const Real _phi_M; // Electrode (metal) potential [V]
};
