
#pragma once

// Include the base class so it can be extended
#include "ADIntegratedBC.h"

// Forward declare the class being created and the validParams function
class miecSEBV;

//** Tafel relation is used for deposition at Li metal/SE interface

class miecSEBV : public ADIntegratedBC
{
public:
  static InputParameters validParams();
  miecSEBV(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;
  
  // The electric potential of electron
  const ADVariableValue & _potEn;

  /// Get parameters from Material system
  const ADMaterialProperty<Real> & _exchange_current;
  const ADMaterialProperty<Real> & _reaction_rate;
  MaterialPropertyName _LiPotElectrode;
  const ADMaterialProperty<Real> & _LiPotEle;
  const Real & _F_RT;

//  usingIntegratedBCMembers;
};
