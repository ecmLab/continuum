#include "ChargeDensity.h"

#include <cstddef>

registerMooseObject("tecm_testApp", ChargeDensity);

InputParameters
ChargeDensity::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Computes total charge density for Phase_2 EDL diagnostics.");

  params.addCoupledVar("concentrations", "Coupled concentration variables contributing to charge density.");
  params.addRequiredParam<std::vector<Real>>("valences", "Valences z_i for each species.");
  params.addRequiredParam<Real>("faraday_constant", "Faraday constant F in C/mol.");
  params.addParam<Real>("fixed_charge_density", 0.0, "Uniform fixed charge density ρ_fixed [C/m³].");
  params.addParam<MaterialPropertyName>(
      "fixed_charge_density_material",
      "Optional material property with spatially varying fixed charge density [C/m³].");

  return params;
}

ChargeDensity::ChargeDensity(const InputParameters & parameters)
  : AuxKernel(parameters),
    _valences(getParam<std::vector<Real>>("valences")),
    _faraday_constant(getParam<Real>("faraday_constant")),
    _fixed_charge_density_material(isParamValid("fixed_charge_density_material")
                                       ? &getMaterialProperty<Real>("fixed_charge_density_material")
                                       : nullptr),
    _fixed_charge_density_constant(getParam<Real>("fixed_charge_density"))
{
  const unsigned int n_conc = coupledComponents("concentrations");
  if (n_conc != _valences.size())
    paramError("valences", "Number of valences must match number of concentration variables.");

  _concentrations.reserve(n_conc);
  for (unsigned int i = 0; i < n_conc; ++i)
    _concentrations.push_back(&coupledValue("concentrations", i));
}

Real
ChargeDensity::computeValue()
{
  Real charge_density = _fixed_charge_density_constant;
  if (_fixed_charge_density_material)
    charge_density += (*_fixed_charge_density_material)[_qp];

  Real mobile_charge = 0.0;
  for (std::size_t i = 0; i < _concentrations.size(); ++i)
    mobile_charge += _valences[i] * (*_concentrations[i])[_qp];

  charge_density += _faraday_constant * mobile_charge;
  return charge_density;
}
