/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   CZMDiffusiveVariableKernelSmallStrain.C
 * Author: srinath
 *
 * Created on January 25, 2022, 8:27 AM
 */

#include "CZMDiffusiveVariableKernelSmallStrain.h"

registerMooseObject("ecmApp", CZMDiffusiveVariableKernelSmallStrain);

InputParameters
CZMDiffusiveVariableKernelSmallStrain::validParams()
{
  InputParameters params = CZMDiffusiveVariableKernelBase::validParams();

  params.addClassDescription(
      "CZM Interface kernel to use when using the Small Strain kinematic formulation.");

  return params;
}

CZMDiffusiveVariableKernelSmallStrain::CZMDiffusiveVariableKernelSmallStrain(const InputParameters & parameters)
  : CZMDiffusiveVariableKernelBase(parameters)
{
}

Real
CZMDiffusiveVariableKernelSmallStrain::computeDResidualDVariable(
    const Moose::DGJacobianType & type) const
{
  Real jac = _dflux_dvariablejump_global[_qp];
  switch (type)
  {
    case Moose::ElementElement: // Residual_sign -1  ddeltaU_ddisp sign -1;
      jac *= _test[_i][_qp] * _phi[_j][_qp];
      break;
    case Moose::ElementNeighbor: // Residual_sign -1  ddeltaU_ddisp sign 1;
      jac *= -_test[_i][_qp] * _phi_neighbor[_j][_qp];
      break;
    case Moose::NeighborElement: // Residual_sign 1  ddeltaU_ddisp sign -1;
      jac *= -_test_neighbor[_i][_qp] * _phi[_j][_qp];
      break;
    case Moose::NeighborNeighbor: // Residual_sign 1  ddeltaU_ddisp sign 1;
      jac *= _test_neighbor[_i][_qp] * _phi_neighbor[_j][_qp];
      break;
  }
  return jac;
}

Real
CZMDiffusiveVariableKernelSmallStrain::computeDResidualDDisplacement(
    const unsigned int & component_j, const Moose::DGJacobianType & type) const
{
  Real jac = _dflux_djump_global[_qp](component_j);
  switch (type)
  {
    case Moose::ElementElement: // Residual_sign -1  ddeltaU_ddisp sign -1;
      jac *= _test[_i][_qp] * _disp_vars[component_j]->phiFace()[_j][_qp];
      break;
    case Moose::ElementNeighbor: // Residual_sign -1  ddeltaU_ddisp sign 1;
      jac *= -_test[_i][_qp] * _disp_vars[component_j]->phiFaceNeighbor()[_j][_qp];
      break;
    case Moose::NeighborElement: // Residual_sign 1  ddeltaU_ddisp sign -1;
      jac *= -_test_neighbor[_i][_qp] * _disp_vars[component_j]->phiFace()[_j][_qp];
      break;
    case Moose::NeighborNeighbor: // Residual_sign 1  ddeltaU_ddisp sign 1;
      jac *= _test_neighbor[_i][_qp] * _disp_vars[component_j]->phiFaceNeighbor()[_j][_qp];
      break;
  }
  return jac;
}
