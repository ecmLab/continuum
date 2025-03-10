
#include "ADMassFluxEEL.h"

registerADMooseObject("ecBetaApp", ADMassFluxEEL);

InputParameters
ADMassFluxEEL::validParams()
{
  InputParameters params = ADKernel::validParams();

  params.addRequiredCoupledVar("Voltage", "The variable representing the voltage.");
  params.addRequiredParam<MaterialPropertyName>("mobility", "The mobility coefficient.");
  params.addRequiredParam<MaterialPropertyName>("conductivity", "The conductivity coefficient.");
  //params.addParam<Real>("zIons", 1, "The charge state of the ions, default be positive 1");
  params.addParam<Real>("c0", "The reference concentration");
  params.addParam<Real>("F", 96485 , "The value of F ,in unit C/mole");
  params.addParam<Real>("R", 8.854 , "The value of R ,in unit J/(mole*K)");
  params.addParam<Real>("T", 298 , "The value of T ,in unit K");
  params.addParam<Real>("scale", 1, "The Scalar Constant, when T = 300K.");
  return params;
}

ADMassFluxEEL::ADMassFluxEEL(const InputParameters & parameters)
   : ADKernel(parameters),
    _grad_V(adCoupledGradient("Voltage")),
    _mobility(getMaterialProperty<Real>("mobility")),
    _conductivity(getMaterialProperty<Real>("conductivity")),
    _F(getParam<Real>("F")),
    _R(getParam<Real>("R")),
    _T(getParam<Real>("T")),
    _c0(getParam<Real>("c0")),
    //_zIons(getParam<Real>("zIons")),
    _scale(getParam<Real>("scale"))


{
}

ADReal
ADMassFluxEEL::computeQpResidual()
{
  //RealVectorValue ion_velocity = _F_RT * _diffusivity[_qp] * _grad_V[_qp];
  ADRealVectorValue MassFluxdDueToConcGrad= _mobility[_qp] * _R * _T * _c0 * _grad_u[_qp] / _test[_i][_qp];
  ADRealVectorValue MassFluxdDueToPotGrad = - _conductivity[_qp] / _F * _grad_V[_qp];
  return _scale * (MassFluxdDueToConcGrad + MassFluxdDueToPotGrad) * _grad_test[_i][_qp];
}
