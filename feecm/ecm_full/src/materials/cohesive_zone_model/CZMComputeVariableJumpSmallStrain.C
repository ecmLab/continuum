/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CZMComputeVariableJumpSmallStrain.C
 * Author: srinath
 * 
 * Created on January 25, 2022, 1:49 PM
 */

#include "CZMComputeVariableJumpSmallStrain.h"


registerMooseObject("TensorMechanicsApp", CZMComputeVariableJumpSmallStrain);

InputParameters
CZMComputeVariableJumpSmallStrain::validParams()
{
  InputParameters params = CZMComputeVariableJumpBase::validParams();
  params.addClassDescription("Compute the total displacement jump across a czm interface in local "
                             "coordinates for the Small Strain kinematic formulation");

  return params;
}

CZMComputeVariableJumpSmallStrain::CZMComputeVariableJumpSmallStrain(
    const InputParameters & parameters)
  : CZMComputeVariableJumpBase(parameters)
{
}

void
CZMComputeVariableJumpSmallStrain::computeLocalVariableJump()
{
  /// This is a scalar variable so there are no rotations to perform
  _interface_variable_jump[_qp] = _variable_jump_global[_qp];
}
