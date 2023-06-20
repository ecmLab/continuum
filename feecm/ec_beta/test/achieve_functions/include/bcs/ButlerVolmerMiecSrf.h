
#pragma once

// Include the base class so it can be extended
#include "ADIntegratedBC.h"

// Forward declare the class being created and the validParams function
class ButlerVolmerMiecSrf : public ADIntegratedBC
{
public:
  static InputParameters validParams();
  ButlerVolmerMiecSrf(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // Get parameters from input
 const Real & _F_RT;
 const Real & _reaction_rate;
 const Real & _ele_conc;
 const Real & _exchange_current;

  // The electric potential
  const ADVariableValue & _potEn;

//  usingIntegratedBCMembers;
};
