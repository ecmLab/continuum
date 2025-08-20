
#include "ScaledCoupledVarNeumannBC.h"

registerMooseObject("ecmApp", ScaledCoupledVarNeumannBC);

// defineLegacyParams(ScaledCoupledVarNeumannBC);

InputParameters
ScaledCoupledVarNeumannBC::validParams()

{
  InputParameters params = IntegratedBC::validParams();
  params.addRequiredCoupledVar("v", "Coupled variable setting the gradient on the boundary.");
  params.addParam<Real>("scale", "Scaling factor for variable to impose BC");
  params.addClassDescription("Imposes the integrated boundary condition "
                             "$\\frac{\\partial u}{\\partial n}=v*scale$, "
                             "where $v$ is a variable.");
  return params;
}

ScaledCoupledVarNeumannBC::ScaledCoupledVarNeumannBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
        _scale(isParamValid("scale") ? getParam<Real>("scale") : 1.0),
        _coupled_var(coupledValue("v"))
{
}

Real
ScaledCoupledVarNeumannBC::computeQpResidual()
{
  return -_test[_i][_qp] * _coupled_var[_qp] * _scale;
}
