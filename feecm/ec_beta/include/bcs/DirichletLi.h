
#pragma once

#include "ADNodalBC.h"

/// Implements a coupled Dirichlet BC where u = alpha * some_var on the boundary.
class DirichletLi : public ADNodalBC
{
public:
  static InputParameters validParams();
  DirichletLi(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

private:
  // The electric potential of Electron
  const ADVariableValue & _potEn;

//  usingNodalBCMembers;
};

