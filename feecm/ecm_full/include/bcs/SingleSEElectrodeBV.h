
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
  const ADMaterialProperty<Real> & _ionic_conductivity;
  const ADMaterialProperty<Real> & _exchange_current;
  const ADMaterialProperty<Real> & _reaction_rate;
  MaterialPropertyName _LiPotElectrode;
  const ADMaterialProperty<Real> & _LiPotEle;
  const Real & _F_RT;

//  usingIntegratedBCMembers;
};
