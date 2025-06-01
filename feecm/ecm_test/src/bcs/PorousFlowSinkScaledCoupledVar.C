
#include "PorousFlowSinkScaledCoupledVar.h"

registerMooseObject("ecmApp", PorousFlowSinkScaledCoupledVar);

InputParameters
PorousFlowSinkScaledCoupledVar::validParams()
{
    InputParameters params = PorousFlowSink::validParams();
    params.addRequiredCoupledVar("v", "Coupled variable setting the gradient on the boundary.");
    params.addParam<Real>("scale", "Scaling factor for variable to impose BC");
    params.addClassDescription("Imposes the integrated boundary condition "
                             "$\\frac{\\partial u}{\\partial n}=v$, "
                             "where $v$ is a variable.");
    return params;
}

PorousFlowSinkScaledCoupledVar::PorousFlowSinkScaledCoupledVar(const InputParameters& parameters)
    : PorousFlowSink(parameters),
    _scale(isParamValid("scale") ? getParam<Real>("scale") : 1.0),
    _coupled_var(coupledValue("v"))

{
}

Real
PorousFlowSinkScaledCoupledVar::computeQpResidual()
{
    return PorousFlowSink::computeQpResidual() * _coupled_var[_qp] * _scale;
}
