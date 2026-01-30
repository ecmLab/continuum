#pragma once

#include "Material.h"

/**
 * SpeciesFluxNoEN implements species flux calculation for TECM without electroneutrality
 * Based on the derivation in TECM_noEN_transference.pdf
 * 
 * Species flux formula:
 * j_i = -∑_j M°_ij ∇μ_j + (t_i)/(z_i F) i
 * 
 * Where:
 * - M°_ij is the mobility matrix (typically diagonal)
 * - t_i is the transference number of species i
 * - z_i is the valence of species i
 * - i is the total current density
 */
class SpeciesFluxNoEN : public Material
{
public:
  static InputParameters validParams();

  SpeciesFluxNoEN(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Names of species flux properties
  const std::vector<MaterialPropertyName> _flux_names;

  /// Concentration variables for each species
  std::vector<const ADVariableValue *> _concentrations;

  /// Chemical potential gradients for each species
  std::vector<const ADVariableGradient *> _grad_mu;

  /// Current density vector
  const ADMaterialProperty<ADRealVectorValue> & _current_density;

  /// Diffusivities for each species
  std::vector<const ADMaterialProperty<Real> *> _diffusivities;

  /// Transference numbers for each species
  std::vector<const ADMaterialProperty<Real> *> _transference_numbers;

  /// Valences for each species
  const std::vector<Real> _valences;

  /// Temperature
  const ADVariableValue & _T;

  /// Physical constants
  const Real _R;  // Gas constant
  const Real _F;  // Faraday constant

  /// Species flux vectors (output)
  std::vector<ADMaterialProperty<ADRealVectorValue> *> _species_fluxes;

  /// Number of species
  const unsigned int _n_species;
};