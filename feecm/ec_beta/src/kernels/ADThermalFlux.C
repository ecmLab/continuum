#include "ADThermalFlux.h"

registerADMooseObject("ecBetaApp", ADThermalFlux);

InputParameters
ADThermalFlux::validParams()
{
    InputParameters params = ADKernel::validParams();
    params.addClassDescription("Compute The Thermal Flux");
    params.addParam<MaterialPropertyName>("k","The Thermal Conductivity");
    params.addParam<Real>("scale",1,"Scaling Parameter");
    return params;
}
ADThermalFlux::ADThermalFlux(const InputParameters & parameters) : ADKernel(parameters),
//
_k(getMaterialProperty<Real>("k")),
_scale(getParam<Real>("scale"))
{
}
ADReal
ADThermalFlux::computeQpResidual()
{
  return _k[_qp]*_scale*_grad_test[_i][_qp] * _grad_u[_qp];

}
