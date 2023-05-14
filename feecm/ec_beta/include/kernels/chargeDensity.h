
#pragma once

#include "ADKernel.h"

class chargeDensity;

class chargeDensity : public ADKernel
{
public:
  static InputParameters validParams();
  chargeDensity(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;
//  virtual Real computeQpJacobian() override;

  const ADVariableValue & _conIons;
//  const VariableGradient & _grad_con;
  const Real & _zIons;
  const Real & _scale;

//  usingGenericKernelMembers;
};
