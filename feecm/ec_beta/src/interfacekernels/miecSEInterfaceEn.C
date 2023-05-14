#include "miecSEInterfaceEn.h"

registerMooseObject("ecBetaApp", miecSEInterfaceEn);

InputParameters
miecSEInterfaceEn::validParams()
{
  InputParameters params = ADInterfaceKernel::validParams();
  // Add a coupled parameter: potEn
    params.addRequiredCoupledVar("potLi", "The potential of Lithium");

// Add a parameter with a default value; this value can be overridden in the input file.
    params.addParam<Real>(
        "F_RT",
        38.68,
        "The constant of F/RT,in unit 1/V, when T = 300K.");

  return params;
}

miecSEInterfaceEn::miecSEInterfaceEn(const InputParameters & parameters)
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
   _F_RT(getParam<Real>("F_RT"))

{
}

ADReal
miecSEInterfaceEn::computeQpResidual(Moose::DGResidualType type)
{
  ADReal r = 0;
  ADReal k1 = std::exp(_reaction_rate[_qp] * _F_RT * (_potLi[_qp] + _u[_qp]));
  ADReal k2 = std::exp(- (1 - _reaction_rate[_qp]) * _F_RT * (_potLi[_qp] + _u[_qp]));
//  Real k1 = std::exp(_reaction_rate[_qp] * _F_RT * (_neighbor_value[_qp] + _u[_qp]));
//  Real k2 = std::exp(- (1 - _reaction_rate[_qp]) * _F_RT * (_neighbor_value[_qp] + _u[_qp]));
//  Real k1 = _electron_concentration[_qp] * _exchange_current[_qp] * _F_RT * (_neighbor_value[_qp] + _u[_qp]);

  switch (type)
  {
    case Moose::Element:
      r = _test[_i][_qp] * (_electron_concentration[_qp] * _exchange_current[_qp] *(k1 - k2)  - 10000*_metal_conductivity[_qp] * _grad_neighbor_value[_qp] * _normals[_qp]);
//      r = _test[_i][_qp] * ( k1 - 10000*_metal_conductivity[_qp] * _grad_neighbor_value[_qp] * _normals[_qp]);
      break;

    case Moose::Neighbor:
      r = _test_neighbor[_i][_qp] * (_electron_concentration[_qp] * _exchange_current[_qp] * (k1 - k2) + 10000*_electronic_conductivity[_qp] * _grad_u[_qp] * _normals[_qp]);
//      r = _test_neighbor[_i][_qp] * (k1 + 10000*_electronic_conductivity[_qp] * _grad_u[_qp] * _normals[_qp]);
      break;
  }

  return r;
}

/**
ADReal
miecSEInterfaceEn::computeQpJacobian(Moose::DGJacobianType type)
{
  ADReal jac = 0;
  ADReal k1 = std::exp(_reaction_rate[_qp] * _F_RT * (_potLi[_qp] + _u[_qp]));
  ADReal k2 = std::exp(- (1 - _reaction_rate[_qp]) * _F_RT * (_potLi[_qp] + _u[_qp]));
//  Real k1 = std::exp(_reaction_rate[_qp] * _F_RT * (_neighbor_value[_qp] + _u[_qp]));
//  Real k2 = std::exp(- (1 - _reaction_rate[_qp]) * _F_RT * (_neighbor_value[_qp] + _u[_qp]));
//  Real k1 = _electron_concentration[_qp] * _exchange_current[_qp] * _F_RT * (_neighbor_value[_qp] + _u[_qp]);

  switch (type)
  {
    case Moose::ElementElement:
      jac = _test[_i][_qp] * _electron_concentration[_qp] * _exchange_current[_qp] * _F_RT * _phi[_j][_qp] *(_reaction_rate[_qp]*k1 + (1-_reaction_rate[_qp])*k2);
//     jac = _test[_i][_qp] * _electron_concentration[_qp] * _exchange_current[_qp] * _F_RT * _phi[_j][_qp];
      break;

    case Moose::NeighborNeighbor:
//      jac = _test_neighbor[_i][_qp] * _electron_concentration[_qp] * _exchange_current[_qp] * _F_RT * _phi_neighbor[_j][_qp] *(_reaction_rate[_qp]*k1 + (1-_reaction_rate[_qp])*k2);
//      jac = _test_neighbor[_i][_qp] * _electron_concentration[_qp] * _exchange_current[_qp] * _F_RT * _phi_neighbor[_j][_qp];
      break;

    case Moose::NeighborElement:
      jac = _test_neighbor[_i][_qp] * (_electron_concentration[_qp] * _exchange_current[_qp] * _F_RT * _phi[_j][_qp]*(_reaction_rate[_qp]*k1 + (1-_reaction_rate[_qp])*k2) + 10000 * _electronic_conductivity[_qp] * _grad_phi[_j][_qp] * _normals[_qp] );
//      jac = _test_neighbor[_i][_qp] * (_electron_concentration[_qp] * _exchange_current[_qp] * _F_RT * _phi[_j][_qp] + 10000 * _electronic_conductivity[_qp] * _grad_phi[_j][_qp] * _normals[_qp] );
      break;

    case Moose::ElementNeighbor:
      jac = _test[_i][_qp] * (- 10000 * _metal_conductivity[_qp] * _grad_phi_neighbor[_j][_qp] * _normals[_qp] );
//      jac = _test[_i][_qp] * (_electron_concentration[_qp] * _exchange_current[_qp] * _F_RT * _phi_neighbor[_j][_qp]*(_reaction_rate[_qp]*k1 + (1-_reaction_rate[_qp])*k2) - 10000 * _metal_conductivity[_qp] * _grad_phi_neighbor[_j][_qp] * _normals[_qp] );
//      jac = _test[_i][_qp] * (_electron_concentration[_qp] * _exchange_current[_qp] * _F_RT * _phi_neighbor[_j][_qp] - 10000 * _metal_conductivity[_qp] * _grad_phi_neighbor[_j][_qp] * _normals[_qp]);
      break;
  }

  return jac;
}
**/
