
#pragma once

// Include the base class so it can be extended
#include "ADIntegratedBC.h"

// Forward declare the class being created and the validParams function
class miecSEEnBV;

//** Tafel relation is used for deposition at Li metal/SE interface

class miecSEEnBV : public ADIntegratedBC
{
public:
  static InputParameters validParams();
  miecSEEnBV(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // The electric potential in metal of pore
  const ADVariableValue & _potMt;
  const ADVariableGradient & _potMt_gradient;

  // The Li+ potential
  const ADVariableValue & _potLi;

  /// Get parameters from Material system
  const ADMaterialProperty<Real> & _metal_conductivity;
  const ADMaterialProperty<Real> & _exchange_current;
  const ADMaterialProperty<Real> & _reaction_rate;
  const ADMaterialProperty<Real> & _electron_concentration;
  const Real & _F_RT;

//  usingIntegratedBCMembers;
};
