#include "PhaseFieldLaplace.h"

registerADMooseObject("ecBetaApp", PhaseFieldLaplace);

InputParameters
PhaseFieldLaplace::validParams()
{
    InputParameters params = ADKernel::validParams();
    params.addClassDescription("Compute Laplace Phase Field");
    params.addParam<MaterialPropertyName>("k0","k0","the diffsuivity coefficient");
    params.addParam<MaterialPropertyName>("scale",1,"It can be Mobility");
    return params;
}
PhaseFieldLaplace::PhaseFieldLaplace(const InputParameters & parameters) : ADKernel(parameters),
//
_k0(getADMaterialProperty<Real>("k0")),
_scale(getADMaterialProperty<Real>("scale"))
{
}
ADReal
PhaseFieldLaplace::computeQpResidual()
{
    ADReal k = _k0[_qp]*_scale[_qp]*_grad_test[_i][_qp] * _grad_u[_qp];
    return k;

}
