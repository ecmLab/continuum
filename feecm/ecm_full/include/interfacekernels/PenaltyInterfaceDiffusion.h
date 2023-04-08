
#pragma once

#include "ADInterfaceKernel.h"

// Forward Declarations
class PenaltyInterfaceDiffusion;

/**
 * DG kernel for interfacing diffusion between two variables on adjacent blocks
 */
class PenaltyInterfaceDiffusion : public ADInterfaceKernel
{
public:
  static InputParameters validParams();
 
  PenaltyInterfaceDiffusion(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual(Moose::DGResidualType type) override;
//  virtual ADReal computeQpJacobian(Moose::DGJacobianType type) override;

  const Real _penalty;
};
