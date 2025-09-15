#include "SpeciesFluxNoEN.h"

registerADMooseObject("tecm_testApp", SpeciesFluxNoEN);

InputParameters
SpeciesFluxNoEN::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Computes species fluxes for TECM without electroneutrality assumption using "
      "transference numbers and explicit valences. Implements: "
      "j_i = -M_i ∇μ_i + (t_i)/(z_i F) i");

  params.addRequiredParam<std::vector<MaterialPropertyName>>("species_fluxes",
                                                             "Names of species flux properties");
  params.addRequiredParam<std::vector<VariableName>>("concentrations",
                                                     "Concentration variables for each species");
  params.addRequiredParam<std::vector<VariableName>>("chemical_potentials",
                                                     "Chemical potential variables for each species");
  params.addRequiredParam<std::vector<MaterialPropertyName>>("diffusivities",
                                                             "Diffusivities for each species");
  params.addRequiredParam<MaterialPropertyName>("current_density",
                                                "Total current density");
  params.addRequiredParam<std::vector<MaterialPropertyName>>("transference_numbers",
                                                             "Transference numbers for each species");
  params.addRequiredParam<std::vector<Real>>("valences",
                                            "Valences for each species (z_i)");
  params.addRequiredParam<VariableName>("temperature", "Temperature variable");
  params.addRequiredParam<Real>("gas_constant", "Universal gas constant");
  params.addRequiredParam<Real>("faraday_constant", "Faraday's constant");

  return params;
}

SpeciesFluxNoEN::SpeciesFluxNoEN(const InputParameters & parameters)
  : Material(parameters),
    _flux_names(getParam<std::vector<MaterialPropertyName>>("species_fluxes")),
    _current_density(getADMaterialProperty<ADRealVectorValue>("current_density")),
    _valences(getParam<std::vector<Real>>("valences")),
    _T(adCoupledValue("temperature")),
    _R(getParam<Real>("gas_constant")),
    _F(getParam<Real>("faraday_constant")),
    _n_species(_flux_names.size())
{
  // Get concentration values
  const auto & c_names = getParam<std::vector<VariableName>>("concentrations");
  _concentrations.resize(_n_species);
  for (unsigned int i = 0; i < _n_species; ++i)
    _concentrations[i] = &adCoupledValue(c_names[i]);

  // Get chemical potential gradients
  const auto & mu_names = getParam<std::vector<VariableName>>("chemical_potentials");
  _grad_mu.resize(_n_species);
  for (unsigned int i = 0; i < _n_species; ++i)
    _grad_mu[i] = &adCoupledGradient(mu_names[i]);

  // Get diffusivities
  const auto & D_names = getParam<std::vector<MaterialPropertyName>>("diffusivities");
  _diffusivities.resize(_n_species);
  for (unsigned int i = 0; i < _n_species; ++i)
    _diffusivities[i] = &getADMaterialProperty<Real>(D_names[i]);

  // Get transference numbers
  const auto & t_names = getParam<std::vector<MaterialPropertyName>>("transference_numbers");
  _transference_numbers.resize(_n_species);
  for (unsigned int i = 0; i < _n_species; ++i)
    _transference_numbers[i] = &getADMaterialProperty<Real>(t_names[i]);

  // Declare species flux properties
  _species_fluxes.resize(_n_species);
  for (unsigned int i = 0; i < _n_species; ++i)
    _species_fluxes[i] = &declareADProperty<ADRealVectorValue>(_flux_names[i]);

  // Validate inputs
  if (c_names.size() != _n_species)
    paramError("concentrations", "Number of concentrations must match number of flux names");
  
  if (mu_names.size() != _n_species)
    paramError("chemical_potentials", "Number of chemical potentials must match number of flux names");
    
  if (D_names.size() != _n_species)
    paramError("diffusivities", "Number of diffusivities must match number of flux names");
    
  if (t_names.size() != _n_species)
    paramError("transference_numbers", "Number of transference numbers must match number of flux names");
    
  if (_valences.size() != _n_species)
    paramError("valences", "Number of valences must match number of flux names");
}

void
SpeciesFluxNoEN::computeQpProperties()
{
  for (unsigned int i = 0; i < _n_species; ++i)
  {
    // Compute mobility: M_i = D_i * c_i / (R * T)
    ADReal mobility = (*_diffusivities[i])[_qp] * std::max((*_concentrations[i])[_qp], 1e-12) / 
                      (_R * _T[_qp]);
    
    // Diffusion contribution: -M_i ∇μ_i
    ADRealVectorValue diffusion_flux = -mobility * (*_grad_mu[i])[_qp];
    
    // Electromigration contribution: (t_i)/(z_i F) i
    ADRealVectorValue migration_flux = ((*_transference_numbers[i])[_qp] / (_valences[i] * _F)) * 
                                       _current_density[_qp];
    
    // Total flux
    (*_species_fluxes[i])[_qp] = diffusion_flux + migration_flux;
  }
}