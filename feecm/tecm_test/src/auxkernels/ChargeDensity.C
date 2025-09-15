#include "ChargeDensity.h"

registerMooseObject("tecm_testApp", ChargeDensity);

InputParameters
ChargeDensity::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Computes charge density for TECM without electroneutrality: "
      "ρₑ = F(z₊c₊ + z₋c₋) + ρfixed");
  
  params.addRequiredCoupledVar("concentrations",
                               "Concentration variables for each species");
  params.addRequiredParam<std::vector<Real>>("valences",
                                            "Valences for each species (z_i)");
  params.addRequiredParam<Real>("faraday_constant", "Faraday's constant");
  params.addParam<Real>("fixed_charge_density", 0.0, "Fixed charge density (default = 0)");

  return params;
}

ChargeDensity::ChargeDensity(const InputParameters & parameters)
  : AuxKernel(parameters),
    _valences(getParam<std::vector<Real>>("valences")),
    _F(getParam<Real>("faraday_constant")),
    _rho_fixed(getParam<Real>("fixed_charge_density")),
    _n_species(coupledComponents("concentrations"))
{
  // Get concentration values
  _concentrations.resize(_n_species);
  for (unsigned int i = 0; i < _n_species; ++i)
    _concentrations[i] = &coupledValue("concentrations", i);

  // Validate inputs
  if (_valences.size() != _n_species)
    paramError("valences", "Number of valences must match number of concentrations");
}

Real
ChargeDensity::computeValue()
{
  // Compute ρₑ = F(z₊c₊ + z₋c₋) + ρfixed
  Real charge_density = _rho_fixed;
  
  for (unsigned int i = 0; i < _n_species; ++i)
  {
    charge_density += _F * _valences[i] * (*_concentrations[i])[_qp];
  }
  
  return charge_density;
}