
#include "PenaltyIntDiffusion.h"

registerMooseObject("ecmApp", PenaltyIntDiffusion);

InputParameters
PenaltyIntDiffusion::validParams()
{
  InputParameters params = ADInterfaceKernel::validParams();
  params.addRequiredParam<Real>(
      "penalty", "The penalty that penalizes jump between master and neighbor variables.");
  return params;
}

PenaltyIntDiffusion::PenaltyIntDiffusion(const InputParameters & parameters)
  : ADInterfaceKernel(parameters), _penalty(getParam<Real>("penalty"))
{
}

ADReal
PenaltyIntDiffusion::computeQpResidual(Moose::DGResidualType type)
{
  ADReal r = 0;

  switch (type)
  {
    case Moose::Element:
      r = _test[_i][_qp] * _penalty * (_u[_qp] - _neighbor_value[_qp]);
      break;

    case Moose::Neighbor:
      r = _test_neighbor[_i][_qp] * -_penalty * (_u[_qp] - _neighbor_value[_qp]);
      break;
  }

  return r;
}
