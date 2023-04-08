
#pragma once

#include "Kernel.h"

class sourceTerm;

class sourceTerm : public Kernel
{
public:
  static InputParameters validParams();
  sourceTerm(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  const MaterialProperty<Real> & _srcCoef;
  const Function & _function;
  const Real & _F_RT;

//  usingGenericKernelMembers;
};
