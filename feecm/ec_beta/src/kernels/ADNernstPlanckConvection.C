
#include "ADNernstPlanckConvection.h"

registerADMooseObject("ecBetaApp", ADNernstPlanckConvection);

InputParameters
ADNernstPlanckConvection::validParams()
{
  InputParameters params = ADKernel::validParams();

  params.addRequiredCoupledVar("Voltage", "The variable representing the voltage.");
  params.addRequiredParam<MaterialPropertyName>("diffusivity", "The diffusivity coefficient.");
  params.addParam<Real>("zIons", 1, "The charge state of the ions, default be positive 1");
  params.addParam<Real>("F_RT", 38.68, "The constant of F/RT,in unit 1/V, when T = 300K.");
  params.addParam<Real>("scale", 1, "The constant of F/RT,in unit 1/V, when T = 300K.");
  return params;
}

ADNernstPlanckConvection::ADNernstPlanckConvection(const InputParameters & parameters)
   : ADKernel(parameters),
    _grad_V(adCoupledGradient("Voltage")),
    _diffusivity(getMaterialProperty<Real>("diffusivity")),
    _scale(getParam<Real>("scale")),
    _F_RT(getParam<Real>("F_RT")),
    _zIons(getParam<Real>("zIons"))

{
}

ADReal
ADNernstPlanckConvection::computeQpResidual()
{
  //RealVectorValue ion_velocity = _F_RT * _diffusivity[_qp] * _grad_V[_qp];
  return  _diffusivity[_qp]*_grad_test[_i][_qp]* _grad_u[_qp] + _diffusivity[_qp]*_zIons * _F_RT * _grad_V[_qp] * _test[_i][_qp]* _grad_u[_qp];
}
