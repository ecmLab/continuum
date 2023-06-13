
#pragma once

#include "AuxKernel.h"

// Forward Declarations
class CurrentDensity;

/**
 * Auxiliary kernel responsible for computing the Current density given
 * several fluid properties and the potential gradient.
 */
class CurrentDensity : public AuxKernel
{
public:
  static InputParameters validParams();

  CurrentDensity(const InputParameters & parameters);

protected:
  /**
   * AuxKernels MUST override computeValue.  computeValue() is called on
   * every quadrature point.  For Nodal Auxiliary variables those quadrature
   * points coincide with the nodes.
   */
  virtual Real computeValue() override;
//  virtual Real computeValue();

  /// Will hold 0, 1, or 2 corresponding to x, y, or z.
  int _component;

  /// The gradient of a coupled variable
  const VariableGradient & _potential_gradient;

  /// Get the conductivity from Material properties
  MaterialPropertyName _conductivity;
  const ADMaterialProperty<Real> & _conductivity_coef;
};
