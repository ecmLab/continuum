#include "ADCurrentFluxEEL.h"

registerADMooseObject("ecBetaApp", ADCurrentFluxEEL);

InputParameters
ADCurrentFluxEEL::validParams()
{
    InputParameters params = ADKernel::validParams();
    params.addClassDescription("Compute the the transport of ions due to migration and diffusion EEL app");
    params.addRequiredCoupledVar("c", "The Concentration field");
    params.addRequiredParam<MaterialPropertyName>("conductivity", "The conductivity coefficient.");
    params.addParam<Real>("F", 964523, "The Value of F");
    params.addParam<Real>("R", 8.845, "The Value of R");
    params.addParam<Real>("T", 298, "The Value of T=300K");
    params.addParam<Real>("c0", "The reference concentration");
    params.addParam<Real>("scale", 1, "The constant of F/RT,in unit 1/V, when T = 300K.");
    return params;
}

ADCurrentFluxEEL::ADCurrentFluxEEL(const InputParameters & parameters)
  : ADKernel(parameters),
  _c(adCoupledValue("c")),
  _grad_c(adCoupledGradient("c")),
  _conductivity(getMaterialProperty<Real>("conductivity")),
   _F(getParam<Real>("F")),
   _R(getParam<Real>("R")),
   _T(getParam<Real>("T")),
  _c0(getParam<Real>("c0")),
  _scale(getParam<Real>("scale"))
{
}

ADReal
ADCurrentFluxEEL::computeQpResidual()
{

  ADRealVectorValue gradMu = _R*_T*_c0 / _c[_qp] * _grad_c[_qp];
  ADRealVectorValue gradPhi =  _grad_u[_qp];
  return _scale * _conductivity[_qp] * _grad_test[_i][_qp] * (gradPhi + gradMu/_F);
}