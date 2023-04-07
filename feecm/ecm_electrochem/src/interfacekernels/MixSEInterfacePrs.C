#include "MixSEInterfacePrs.h"

registerMooseObject("ecmElectrochemApp", MixSEInterfacePrs);

InputParameters
MixSEInterfacePrs::validParams()
{
  InputParameters params = ADInterfaceKernel::validParams();
  // Add a coupled parameter: potEn
    params.addRequiredCoupledVar("potLi", "The potential of Lithium");

// Add the required parameter the input file.
    params.addRequiredParam<Real>("pressure", "The pressure developed in the void,in unit MPa.");

// Add a parameter with a default value; this value can be overridden in the input file.
    params.addParam<Real>(
        "F_RT",
        38.68,
        "The constant of F/RT,in unit 1/V, when T = 300K.");

// Add a parameter with a default value; this value can be overridden in the input file.
    params.addParam<Real>(
        "Vm_RT",
        0.0052,
        "The constant of Vm/RT,in unit 1/MPa, when T = 300K.");

  return params;
}

MixSEInterfacePrs::MixSEInterfacePrs(const InputParameters & parameters)
  : ADInterfaceKernel(parameters),

  // Couple to the potential of Li
   _potLi(coupledValue("potLi")),

  // Get the gradient of the variable
   _potLi_gradient(coupledGradient("potLi")),

 // Get the parameters from the input file
   _electron_concentration(getADMaterialProperty<Real>("electron_concentration")),
   _ionic_conductivity(getADMaterialProperty<Real>("ionic_conductivity")),
   _electronic_conductivity(getADMaterialProperty<Real>("electronic_conductivity")),
   _metal_conductivity(getADMaterialProperty<Real>("metal_conductivity")),
   _exchange_current(getADMaterialProperty<Real>("exchange_current")),
   _reaction_rate(getADMaterialProperty<Real>("reaction_rate")),
   _pressure(getParam<Real>("pressure")),
   _F_RT(getParam<Real>("F_RT")),
   _Vm_RT(getParam<Real>("Vm_RT"))

{
}

ADReal
MixSEInterfacePrs::computeQpResidual(Moose::DGResidualType type)
{
  ADReal r = 0;
  ADReal k1 = std::exp(_reaction_rate[_qp] * (_F_RT*(_potLi[_qp] + _u[_qp]) - _Vm_RT*_pressure));
  ADReal k2 = std::exp(- (1 - _reaction_rate[_qp]) * (_F_RT*(_potLi[_qp] + _u[_qp]) - _Vm_RT*_pressure));

  switch (type)
  {
    case Moose::Element:
      r = _test[_i][_qp] * (_electron_concentration[_qp] * _exchange_current[_qp] *(k1 - k2)  - 10000*_metal_conductivity[_qp] * _grad_neighbor_value[_qp] * _normals[_qp]);
      break;

    case Moose::Neighbor:
      r = _test_neighbor[_i][_qp] * (_electron_concentration[_qp] * _exchange_current[_qp] * (k1 - k2) + 10000*_electronic_conductivity[_qp] * _grad_u[_qp] * _normals[_qp]);
      break;
  }

  return r;
}

/**
ADReal
MixSEInterfacePrs::computeQpJacobian(Moose::DGJacobianType type)
{
  ADReal jac = 0;
  ADReal k1 = std::exp(_reaction_rate[_qp] * (_F_RT*(_potLi[_qp] + _u[_qp]) - _Vm_RT*_pressure));
  ADReal k2 = std::exp(- (1 - _reaction_rate[_qp]) * (_F_RT*(_potLi[_qp] + _u[_qp]) - _Vm_RT*_pressure));

  switch (type)
  {
    case Moose::ElementElement:
      jac = _test[_i][_qp] * _electron_concentration[_qp] * _exchange_current[_qp] * _F_RT * _phi[_j][_qp] *(_reaction_rate[_qp]*k1 + (1-_reaction_rate[_qp])*k2);
      break;

    case Moose::NeighborNeighbor:
      break;

    case Moose::NeighborElement:
      jac = _test_neighbor[_i][_qp] * (_electron_concentration[_qp] * _exchange_current[_qp] * _F_RT * _phi[_j][_qp]*(_reaction_rate[_qp]*k1 + (1-_reaction_rate[_qp])*k2) + 10000 * _electronic_conductivity[_qp] * _grad_phi[_j][_qp] * _normals[_qp] );
      break;

    case Moose::ElementNeighbor:
      jac = _test[_i][_qp] * (- 10000 * _metal_conductivity[_qp] * _grad_phi_neighbor[_j][_qp] * _normals[_qp] );
      break;
  }

  return jac;
}
**/
