
#include "ParamDiffusion.h"

registerADMooseObject("ecmElectrochemApp", ParamDiffusion);

InputParameters
ParamDiffusion::validParams()
{
    InputParameters params = ADKernel::validParams();
    params.addClassDescription("Compute the ionic/electronic conduction.");
    params.addRequiredParam<MaterialPropertyName>("conductivity", "The conductivity, in mS/cm.");
    return params;
}

ParamDiffusion::ParamDiffusion(const InputParameters & parameters)
  : ADKernel(parameters),
    _conductivity(parameters.get<MaterialPropertyName>("conductivity")),
    _conductivity_coef(getADMaterialProperty<Real>(_conductivity))
{
}

ADReal
ParamDiffusion::computeQpResidual()
{
//  return 10000 * _conductivity_coef[_qp] * ADDiffusion<compute_stage>::precomputeQpResidual();
  return 10000 * _conductivity_coef[_qp] * _grad_test[_i][_qp] * _grad_u[_qp];
}

