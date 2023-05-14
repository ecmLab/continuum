
#pragma once

// Include the base class so it can be extended
#include "ADIntegratedBC.h"

// Forward declare the class being created and the validParams function
class SingleSEElectrodeNeumann;

/**
 * An IntegratedBC representing the Neumann boundary condition for mixed conduction.
 */
class SingleSEElectrodeNeumann : public ADIntegratedBC
{
public:
  static InputParameters validParams();
  SingleSEElectrodeNeumann(const InputParameters & parameters);

protected:
  /**
   * This is called to integrate the residual across the boundary.
   */
  virtual ADReal computeQpResidual() override;
 
  /// Get parameters from Material system
  const ADMaterialProperty<Real> & _applied_current;

//  usingIntegratedBCMembers;
};
