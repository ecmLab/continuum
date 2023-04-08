
#include "ParamDiffusion.h"

registerADMooseObject("liExpulsionApp", ParamDiffusion);

InputParameters
ParamDiffusion::validParams()
{
    InputParameters params = ADKernel::validParams();
    params.addClassDescription("Compute the laplacian equation.");
    params.addRequiredParam<MaterialPropertyName>("diffusivity", "The diffusivity coefficient.");
    params.addParam<Real>("scale", 1, "The scale of the equation");
    return params;
}

ParamDiffusion::ParamDiffusion(const InputParameters & parameters)
  : ADKernel(parameters),
    _diffusivity(getMaterialProperty<Real>("diffusivity")),
    _scale(getParam<Real>("scale"))
{
}

ADReal
ParamDiffusion::computeQpResidual()
{
  return _scale * _diffusivity[_qp] * _grad_test[_i][_qp] * _grad_u[_qp];
}

