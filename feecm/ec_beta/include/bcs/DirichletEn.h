
#pragma once

#include "ADNodalBC.h"
class DirichletEn;

/// Implements a coupled Dirichlet BC where u = alpha * some_var on the boundary.
class DirichletEn : public ADNodalBC
{
public:
  static InputParameters validParams();

  DirichletEn(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

private:
  // The electric potential of Li-ions
  const ADVariableValue & _potLi;

//  usingNodalBCMembers;
};

