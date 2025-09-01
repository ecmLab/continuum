
#include "DarcyPressure.h"

registerMooseObject("mece789App", DarcyPressure);

InputParameters
DarcyPressure::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Compute the diffusion term for Darcy pressure ($p$) equation: "
                             "$-\\nabla \\cdot \\frac{\\mathbf{K}}{\\mu} \\nabla p = 0$");
  return params;
}

DarcyPressure::DarcyPressure(const InputParameters & parameters)
  : ADKernel(parameters),
    _permeability(getADMaterialProperty<Real>("permeability")),
    _viscosity(getADMaterialProperty<Real>("viscosity"))
{
}

ADReal
DarcyPressure::computeQpResidual()
{
  return (_permeability[_qp] / _viscosity[_qp]) * _grad_test[_i][_qp] * _grad_u[_qp];
}
