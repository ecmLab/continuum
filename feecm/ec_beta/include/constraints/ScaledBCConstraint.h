#pragma once

#include "ADMortarConstraint.h"


class ScaledBCConstraint;

/**
* Scales an existing contact problem variable to apply on another variable
* The lagrange multiplier of the existing contact problem is used to compute 
* scaled integrated BC on the opposite side of the interface. 
* This assumes that the primary side is always the one on which the bc should 
* be applied on and the secondary side already contains the lagrange_multiplier
*/

class ScaledBCConstraint : public ADMortarConstraint
{
public:

  static InputParameters validParams();
  ScaledBCConstraint(const InputParameters & parameters);

protected:

  virtual ADReal computeQpResidual(Moose::MortarType mortar_type) override;
  
  ADReal _scale;
  bool _primary;

};
