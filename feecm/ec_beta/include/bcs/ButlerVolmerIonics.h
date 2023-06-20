#pragma once

// Include the base class so it can be extended
#include "ADIntegratedBC.h"

// Forward declare the class being created and the validParams function

class ButlerVolmerIonics : public ADIntegratedBC
{
public:
  static InputParameters validParams();
  ButlerVolmerIonics(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;
  
  /// Get parameters from Material system
  const Real & _F_RT;
  const Real & _reaction_rate;
  const Real & _LiPotRef;
  const Real & _exchange_current;
//  MaterialPropertyName _ex_current;
//  const ADMaterialProperty<Real> & _exchange_current;

//  usingIntegratedBCMembers;
};
