
#pragma once

// Include the base class so it can be extended
#include "ADIntegratedBC.h"

// Forward declare the class being created and the validParams function
class miecSEVoidBV;

//** Tafel relation is used for deposition at Li metal/SE interface

class miecSEVoidBV : public ADIntegratedBC
{
public:
  static InputParameters validParams();
  miecSEVoidBV(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // The electric potential
  const ADVariableValue & _potRef;

  /// Get parameters from Material system
  const ADMaterialProperty<Real> & _exchange_current;
  const ADMaterialProperty<Real> & _reaction_rate;
  const ADMaterialProperty<Real> & _electron_concentration;
  const Real & _F_RT;

//  usingIntegratedBCMembers;
};
