#include "ADPoissonChargeRHS.h"

registerMooseObject("tecm_testApp", ADPoissonChargeRHS);

InputParameters
ADPoissonChargeRHS::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addRequiredCoupledVar("c_plus",  "Cation concentration");
  params.addRequiredCoupledVar("c_minus", "Anion concentration");
  params.addParam<Real>("F", 96485.3329,  "Faraday constant [C/mol]");
  params.addParam<Real>("z_plus",  1.0,   "Valence of cation");
  params.addParam<Real>("z_minus", -1.0,  "Valence of anion");
  return params;
}

ADPoissonChargeRHS::ADPoissonChargeRHS(const InputParameters & params)
  : ADKernel(params),
    _c_plus(adCoupledValue("c_plus")),
    _c_minus(adCoupledValue("c_minus")),
    _F(getParam<Real>("F")),
    _z_plus(getParam<Real>("z_plus")),
    _z_minus(getParam<Real>("z_minus"))
{
}

ADReal
ADPoissonChargeRHS::computeQpResidual()
{
  ADReal rho_e = _F * (_z_plus * _c_plus[_qp] + _z_minus * _c_minus[_qp]);
  // Weak form term: - ∫ w * rho_e
  return - _test[_i][_qp] * rho_e;
}
