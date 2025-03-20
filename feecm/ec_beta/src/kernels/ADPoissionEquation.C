#include "ADPoissonEquation.h"

registerMooseObject("ecBetaApp", ADPoissonEquation);

InputParameters
ADPoissonEquation::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Implements the Poisson equation ∇ · (σ/F ∇φ) = zc");
  params.addRequiredCoupledVar("coupled_concentration", "The coupled concentration variable");
  params.addParam<Real>("F", 96485.3383, "Faraday constant in C/mol");
  params.addRequiredParam<Real>("sigma", "Conductivity");
  params.addRequiredParam<Real>("z", "Charge number");
  params.addParam<Real>("scale",1, "Scale Parameter");
  return params;
}

ADPoissonEquation::ADPoissonEquation(const InputParameters & parameters)
  : ADKernel(parameters),
    _coupled_concentration(adCoupledValue("coupled_concentration")),
    _F(getParam<Real>("F")),
    _sigma(getParam<Real>("sigma")),
    _z(getParam<Real>("z")),
    _scale(getParam<Real>("scale"))
{
}

ADReal
ADPoissonEquation::computeQpResidual()
{
  return _sigma / _F *_scale* _grad_test[_i][_qp] * _grad_u[_qp] -
         _test[_i][_qp] * _z * _coupled_concentration[_qp];
}
