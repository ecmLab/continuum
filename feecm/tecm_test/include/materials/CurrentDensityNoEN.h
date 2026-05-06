#pragma once

#include "Material.h"

/**
 * CurrentDensityNoEN implements current density calculation for TECM without electroneutrality
 * Based on the derivation in TECM_noEN_transference.pdf
 * 
 * Current partition formula:
 * i = -σ°∇Φ - (σ°/F)(t₊∇μ₊ + t₋∇μ₋)
 * 
 * Where:
 * - σ° is the reference conductivity
 * - t₊, t₋ are transference numbers (t₊ + t₋ = 1)
 * - μ₊, μ₋ are chemical potentials of cation and anion
 */
class CurrentDensityNoEN : public Material
{
public:
  static InputParameters validParams();

  CurrentDensityNoEN(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Name of the current density property
  const MaterialPropertyName _current_density_name;

  /// Electric potential variable
  const ADVariableGradient & _grad_Phi;

  /// Chemical potential gradients for each species
  std::vector<const ADVariableGradient *> _grad_mu;

  /// Reference electrical conductivity
  const ADMaterialProperty<Real> & _sigma_ref;

  /// Transference numbers for each species
  std::vector<const ADMaterialProperty<Real> *> _transference_numbers;

  /// Valences for each species
  const std::vector<Real> _valences;

  /// Faraday constant
  const Real _F;

  /// Current density vector (output)
  ADMaterialProperty<ADRealVectorValue> & _current_density;

  /// Derivatives
  ADMaterialProperty<ADRealVectorValue> & _d_current_d_grad_Phi;
  std::vector<ADMaterialProperty<ADRealVectorValue> *> _d_current_d_grad_mu;

  /// Number of species
  const unsigned int _n_species;
};