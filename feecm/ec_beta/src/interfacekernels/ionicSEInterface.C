#include "IonicSEInterface.h"

registerMooseObject("ecBetaApp", IonicSEInterface);

InputParameters
IonicSEInterface::validParams()
{
  InputParameters params = ADInterfaceKernel::validParams();

// Add a parameter with a default value; this value can be overridden in the input file.
    params.addParam<Real>(
        "F_RT",
        38.68,
        "The constant of F/RT,in unit V, when T = 300K.");

  return params;
}

IonicSEInterface::IonicSEInterface(const InputParameters & parameters)
  : ADInterfaceKernel(parameters),

 // Get the parameters from the input file
   _ionic_conductivity(getADMaterialProperty<Real>("ionic_conductivity")),
   _metal_conductivity(getADMaterialProperty<Real>("metal_conductivity")),
   _exchange_current(getADMaterialProperty<Real>("exchange_current")),
   _reaction_rate(getADMaterialProperty<Real>("reaction_rate")),
   _F_RT(getParam<Real>("F_RT"))

{
}

ADReal
IonicSEInterface::computeQpResidual(Moose::DGResidualType type)
{
  ADReal r = 0;
  ADReal k1 = std::exp(_reaction_rate[_qp] * _F_RT * (_u[_qp] +  _neighbor_value[_qp]));
  ADReal k2 = std::exp(- (1 - _reaction_rate[_qp]) * _F_RT * (_u[_qp] + _neighbor_value[_qp]));

  switch (type)
  {
    case Moose::Element:
      r = _test[_i][_qp] * 10000*_metal_conductivity[_qp] * _grad_neighbor_value[_qp] * _normals[_qp];
//      r = _test[_i][_qp] * _exchange_current[_qp] * (k1 - k2);
      break;

    case Moose::Neighbor:
//      r = _test_neighbor[_i][_qp] * 10000*_ionic_conductivity[_qp] * _grad_u[_qp] * _normals[_qp];
      r = _test_neighbor[_i][_qp] * _exchange_current[_qp] * (k1 - k2);
      break;
  }

  return r;
}

