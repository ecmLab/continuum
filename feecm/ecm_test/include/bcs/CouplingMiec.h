
#pragma once

// Include the base class so it can be extended
#include "ADIntegratedBC.h"

class CouplingMiec : public ADIntegratedBC
{
public:
  static InputParameters validParams();
  CouplingMiec(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // Get parameters from input
  const Real & _liCrtRef;

  // Get parameters from Material system
  const ADMaterialProperty<Real> & _ionic_conductivity;

  // The Li+ potential
  const ADVariableValue & _potLi;
  const ADVariableGradient & _potLi_gradient;


//  usingIntegratedBCMembers;
};
