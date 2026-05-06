#pragma once

#include "AuxKernel.h"

#include <vector>

/**
 * Computes the mobile + fixed charge density ρ_e = F Σ z_i c_i + ρ_fixed.
 * Useful for diagnostics when solving the Phase_2 Poisson equation.
 */
class ChargeDensity : public AuxKernel
{
public:
  static InputParameters validParams();

  ChargeDensity(const InputParameters & parameters);

protected:
  Real computeValue() override;

private:
  std::vector<const VariableValue *> _concentrations;
  std::vector<Real> _valences;
  const Real _faraday_constant;
  const MaterialProperty<Real> * _fixed_charge_density_material;
  const Real _fixed_charge_density_constant;
};
