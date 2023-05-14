
#pragma once

// Include the base class so it can be extended
#include "ADIntegratedBC.h"

// Forward declare the class being created and the validParams function
class ionicSELiBV;

//** Tafel relation is used for deposition at Li metal/SE interface

class ionicSELiBV : public ADIntegratedBC
{
public:
  static InputParameters validParams();
  ionicSELiBV(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // The electric potential
  const ADVariableValue & _potLi;

  /// Get parameters from Material system
  const ADMaterialProperty<Real> & _exchange_current;
  const ADMaterialProperty<Real> & _reaction_rate;
  const Real & _F_RT;

//  usingIntegratedBCMembers;
};
