
#pragma once

#include "ADInterfaceKernel.h"

class ADInterfaceDiffusion : public ADInterfaceKernel
{
public:
  static InputParameters validParams();

 ADInterfaceDiffusion(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual(Moose::DGResidualType type) override;
//  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

  const ADMaterialProperty<Real> & _D;
  const ADMaterialProperty<Real> & _D_neighbor;
};
