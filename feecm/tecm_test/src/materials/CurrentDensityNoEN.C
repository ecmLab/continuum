#include "CurrentDensityNoEN.h"

registerADMooseObject("tecm_testApp", CurrentDensityNoEN);

InputParameters
CurrentDensityNoEN::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Computes current density for TECM without electroneutrality assumption using "
      "transference numbers and explicit valences. Implements: "
      "i = -σ°∇Φ - (σ°/F)(t₊∇μ₊ + t₋∇μ₋)");

  params.addRequiredParam<MaterialPropertyName>("current_density", 
                                                "Name of the current density property");
  params.addRequiredParam<VariableName>("electric_potential", 
                                        "The electric potential variable");
  params.addRequiredParam<std::vector<VariableName>>("chemical_potentials",
                                                     "Chemical potential variables for each species");
  params.addRequiredParam<std::vector<MaterialPropertyName>>("transference_numbers",
                                                             "Transference numbers for each species");
  params.addRequiredParam<std::vector<Real>>("valences", 
                                            "Valences for each species (z_i)");
  params.addRequiredParam<MaterialPropertyName>("electric_conductivity",
                                                "Reference electrical conductivity");
  params.addRequiredParam<Real>("faraday_constant", "Faraday's constant");

  return params;
}

CurrentDensityNoEN::CurrentDensityNoEN(const InputParameters & parameters)
  : Material(parameters),
    _current_density_name(getParam<MaterialPropertyName>("current_density")),
    _grad_Phi(adCoupledGradient("electric_potential")),
    _sigma_ref(getADMaterialProperty<Real>("electric_conductivity")),
    _valences(getParam<std::vector<Real>>("valences")),
    _F(getParam<Real>("faraday_constant")),
    _current_density(declareADProperty<ADRealVectorValue>(_current_density_name)),
    _d_current_d_grad_Phi(declareADProperty<ADRealVectorValue>("d_" + _current_density_name + "_d_grad_Phi")),
    _n_species(getParam<std::vector<VariableName>>("chemical_potentials").size())
{
  // Get chemical potential gradients
  const auto & mu_names = getParam<std::vector<VariableName>>("chemical_potentials");
  _grad_mu.resize(_n_species);
  for (unsigned int i = 0; i < _n_species; ++i)
    _grad_mu[i] = &adCoupledGradient(mu_names[i]);

  // Get transference numbers
  const auto & t_names = getParam<std::vector<MaterialPropertyName>>("transference_numbers");
  _transference_numbers.resize(_n_species);
  for (unsigned int i = 0; i < _n_species; ++i)
    _transference_numbers[i] = &getADMaterialProperty<Real>(t_names[i]);

  // Initialize derivatives w.r.t. chemical potential gradients
  _d_current_d_grad_mu.resize(_n_species);
  for (unsigned int i = 0; i < _n_species; ++i)
    _d_current_d_grad_mu[i] = &declareADProperty<ADRealVectorValue>(
        "d_" + _current_density_name + "_d_grad_mu_" + std::to_string(i));

  // Validate inputs
  if (_valences.size() != _n_species)
    paramError("valences", "Number of valences must match number of chemical potentials");
  
  if (t_names.size() != _n_species)
    paramError("transference_numbers", "Number of transference numbers must match number of chemical potentials");
}

void
CurrentDensityNoEN::computeQpProperties()
{
  // Initialize current density with Ohmic contribution: i = -σ°∇Φ
  _current_density[_qp] = -_sigma_ref[_qp] * _grad_Phi[_qp];
  
  // Add chemical potential contributions: i -= (σ°/F)(t₊∇μ₊ + t₋∇μ₋)
  ADRealVectorValue mu_contribution(0.0, 0.0, 0.0);
  for (unsigned int i = 0; i < _n_species; ++i)
  {
    mu_contribution += (*_transference_numbers[i])[_qp] * (*_grad_mu[i])[_qp];
  }
  _current_density[_qp] -= (_sigma_ref[_qp] / _F) * mu_contribution;

  // Compute derivatives
  _d_current_d_grad_Phi[_qp] = -_sigma_ref[_qp] * ADRealVectorValue(1.0, 1.0, 1.0);
  
  for (unsigned int i = 0; i < _n_species; ++i)
  {
    (*_d_current_d_grad_mu[i])[_qp] = -(_sigma_ref[_qp] / _F) * (*_transference_numbers[i])[_qp] * 
                                      ADRealVectorValue(1.0, 1.0, 1.0);
  }
}