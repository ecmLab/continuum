

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
