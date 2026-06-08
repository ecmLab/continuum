
#pragma once

#include "AuxKernel.h"

/**
 * Auxiliary kernel responsible for computing the overpotentials: potLi - potEn
 */
class LiPotMiec : public AuxKernel
{
public:
  static InputParameters validParams();

  LiPotMiec(const InputParameters & parameters);

protected:
  /**
   * AuxKernels MUST override computeValue.  computeValue() is called on
   * every quadrature point.  For Nodal Auxiliary variables those quadrature
   * points coincide with the nodes.
   */
  virtual Real computeValue() override;

 // The electric potential of Li-ions
  const VariableValue & _potLi;

 // The electronic potential of Electron
  const VariableValue & _potEn;

};
