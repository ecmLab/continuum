#include "EDLProperties.h"

registerADMooseObject("tecm_testApp", EDLProperties);

InputParameters EDLProperties::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addRequiredParam<Real>("D_plus", "Cation diffusivity D_plus [m^2/s]");
  params.addRequiredParam<Real>("D_minus", "Anion diffusivity D_minus [m^2/s]");
  params.addRequiredParam<Real>("permittivity", "Permittivity ε [F/m]");
  params.addClassDescription("Constant material properties for diffusion coefficients and permittivity in the Phase 2 EDL model.");
  return params;
}

EDLProperties::EDLProperties(const InputParameters & parameters)
  : ADMaterial(parameters),
    _D_plus_value(getParam<Real>("D_plus")),
    _D_minus_value(getParam<Real>("D_minus")),
    _permittivity_value(getParam<Real>("permittivity")),
    _D_plus(declareADProperty<Real>("D_plus")),
    _D_minus(declareADProperty<Real>("D_minus")),
    _permittivity(declareADProperty<Real>("permittivity"))
{
}

void EDLProperties::computeQpProperties()
{
  _D_plus[_qp] = _D_plus_value;
  _D_minus[_qp] = _D_minus_value;
  _permittivity[_qp] = _permittivity_value;
}

