
#include "ADInterfaceDiffusion.h"

registerMooseObject("ecmApp", ADInterfaceDiffusion);

InputParameters
ADInterfaceDiffusion::validParams()
{
  InputParameters params = ADInterfaceKernel::validParams();
  params.addParam<MaterialPropertyName>("D", "D", "The diffusion coefficient.");
  params.addParam<MaterialPropertyName>(
      "D_neighbor", "D_neighbor", "The neighboring diffusion coefficient.");
  params.addClassDescription(
      "The kernel is utilized to establish flux equivalence on an interface for variables.");
  return params;
}

ADInterfaceDiffusion::ADInterfaceDiffusion(const InputParameters & parameters)
  : ADInterfaceKernel(parameters),
    _D(getADMaterialProperty<Real>("D")),
    _D_neighbor(getNeighborADMaterialProperty<Real>("D_neighbor"))
{
}

ADReal
ADInterfaceDiffusion::computeQpResidual(Moose::DGResidualType type)
{
  ADReal r = 0;

  switch (type)
  {
    case Moose::Element:
      r = _test[_i][_qp] * -_D_neighbor[_qp] * _grad_neighbor_value[_qp] * _normals[_qp];
      break;

    case Moose::Neighbor:
      r = _test_neighbor[_i][_qp] * _D[_qp] * _grad_u[_qp] * _normals[_qp];
      break;
  }

  return r;
}
