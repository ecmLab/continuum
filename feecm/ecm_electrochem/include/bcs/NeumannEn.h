
#pragma once

// Include the base class so it can be extended
#include "ADIntegratedBC.h"

// Forward declare the class being created and the validParams function
class NeumannEn;

/**
 * An IntegratedBC representing the Neumann boundary condition for mixed conduction.
 */
class NeumannEn : public ADIntegratedBC
{
public:
  static InputParameters validParams();
  NeumannEn(const InputParameters & parameters);

protected:
  /**
   * This is called to integrate the residual across the boundary.
   */
  virtual ADReal computeQpResidual() override;
 
  // The electric potential of Li-ions
  const ADVariableValue & _potLi;
  const ADVariableGradient & _potLi_grad;

  /// Get parameters from Material system
  const ADMaterialProperty<Real> & _inlet_current;
  const ADMaterialProperty<Real> & _exchange_current;
  const ADMaterialProperty<Real> & _reaction_rate;
  MaterialPropertyName _LiPotElectrode;
  const ADMaterialProperty<Real> & _LiPotEle;
  const Real & _F_RT;

//  usingIntegratedBCMembers;
};
