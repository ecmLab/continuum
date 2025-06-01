
#pragma once

#include "ADInterfaceKernel.h"

// Forward Declarations
class PenaltyIntDiffusion : public ADInterfaceKernel
{
public:
  static InputParameters validParams();

  PenaltyIntDiffusion(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual(Moose::DGResidualType type) override;
//  virtual ADReal computeQpJacobian(Moose::DGJacobianType type) override;

  const Real & _penalty;
};
