
#pragma once

// Include the base class so it can be extended
#include "ADIntegratedBC.h"

// Forward declare the class being created and the validParams function
class NeumannLi;

/**
 * An IntegratedBC representing the Neumann boundary condition for mixed conduction.
 */
class NeumannLi : public ADIntegratedBC
{
public:
  static InputParameters validParams();
  NeumannLi(const InputParameters & parameters);

protected:
  /**
   * This is called to integrate the residual across the boundary.
   */
  virtual ADReal computeQpResidual() override;
 
  /// The gradient of potEn
  const ADVariableGradient & _potEn_grad;

  /// Get parameters from Material system
  const ADMaterialProperty<Real> & _electronic_conductivity;
  const ADMaterialProperty<Real> & _inlet_current;

//  usingIntegratedBCMembers;
};
