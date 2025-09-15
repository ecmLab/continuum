#pragma once

#include "AuxKernel.h"

/**
 * ChargeDensity computes the charge density for TECM-no-EN framework
 * Based on the derivation in TECM_noEN_transference.pdf
 * 
 * Charge density formula:
 * ρₑ = F(z₊c₊ + z₋c₋) + ρfixed
 * 
 * Where:
 * - F is Faraday's constant
 * - z_i are the valences of each species
 * - c_i are the concentrations of each species
 * - ρfixed is any fixed charge density (default = 0)
 */
class ChargeDensity : public AuxKernel
{
public:
  static InputParameters validParams();

  ChargeDensity(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Concentration variables for each species
  std::vector<const VariableValue *> _concentrations;

  /// Valences for each species
  const std::vector<Real> _valences;

  /// Faraday's constant
  const Real _F;

  /// Fixed charge density (optional)
  const Real _rho_fixed;

  /// Number of species
  const unsigned int _n_species;
};