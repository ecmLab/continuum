
#include "PenaltyInterfaceDiffusion.h"

registerMooseObject("electro_chemo_mechApp", PenaltyInterfaceDiffusion);

// InputParameters
// PenaltyInterfaceDiffusion::validParams()
template <>
InputParameters
validParams<PenaltyInterfaceDiffusion>()
{
  InputParameters params = ADInterfaceKernel::validParams();
  params.addRequiredParam<Real>(
      "penalty", "The penalty that penalizes jump between master and neighbor variables.");
  return params;
}

PenaltyInterfaceDiffusion::PenaltyInterfaceDiffusion(const InputParameters & parameters)
  : ADInterfaceKernel(parameters), _penalty(getParam<Real>("penalty"))
{
}

ADReal
PenaltyInterfaceDiffusion::computeQpResidual(Moose::DGResidualType type)
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

/**
ADReal
PenaltyInterfaceDiffusion::computeQpJacobian(Moose::DGJacobianType type)
{
  ADReal jac = 0;

  switch (type)
  {
    case Moose::ElementElement:
      jac = _test[_i][_qp] * _penalty * _phi[_j][_qp];
      break;

    case Moose::ElementNeighbor:
      jac = _test[_i][_qp] * _penalty * -_phi_neighbor[_j][_qp];
      break;

    case Moose::NeighborElement:
      jac = _test_neighbor[_i][_qp] * -_penalty * _phi[_j][_qp];
      break;

    case Moose::NeighborNeighbor:
      jac = _test_neighbor[_i][_qp] * -_penalty * -_phi_neighbor[_j][_qp];
      break;
  }

  return jac;
}
**/
