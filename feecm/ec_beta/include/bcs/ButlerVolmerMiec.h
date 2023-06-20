
#pragma once

// Include the base class so it can be extended
#include "ADIntegratedBC.h"

// Forward declare the class being created and the validParams function
class ButlerVolmerMiec : public ADIntegratedBC
{
public:
  static InputParameters validParams();
  ButlerVolmerMiec(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;
  
  // Get parameters from input
   const Real & _F_RT;
   const Real & _reaction_rate;
   const Real & _LiCrtRef;
   const Real & _LiPotRef;
   const Real & _exchange_current;

   // The electric potential of electron
  const ADVariableValue & _potEn;

//  usingIntegratedBCMembers;
};
