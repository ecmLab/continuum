
#pragma once

// Include the base class so it can be extended
#include "ADIntegratedBC.h"

// Forward declare the class being created and the validParams function
class SingleSEElectrodeBV;

//** Tafel relation is used for deposition at Li metal/SE interface

class SingleSEElectrodeBV : public ADIntegratedBC
{
public:
  static InputParameters validParams();
  SingleSEElectrodeBV(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;
  
  /// Get parameters from Material system
  const Real & _F_RT;
  const Real & _reaction_rate;
  const Real & _LiPotEle;
  const Real & _exchange_current;

//  usingIntegratedBCMembers;
};
