#pragma once

#include "ADInterfaceKernel.h"

// Forward Declarations
class MiecSEInterfacePrs : public ADInterfaceKernel
{
public:
  static InputParameters validParams();

  MiecSEInterfacePrs(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual(Moose::DGResidualType type) override;
//  virtual ADReal computeQpJacobian(Moose::DGJacobianType type) override;
  
    // The electric potential of Li-ions
  const VariableValue & _potLi;
  const VariableGradient & _potLi_gradient;

  /// Get parameters from Material system
  const ADMaterialProperty<Real> & _electron_concentration;
  const ADMaterialProperty<Real> & _ionic_conductivity;
  const ADMaterialProperty<Real> & _electronic_conductivity;
  const ADMaterialProperty<Real> & _metal_conductivity;
  const ADMaterialProperty<Real> & _exchange_current;
  const ADMaterialProperty<Real> & _reaction_rate;
  const Real & _pressure;
  const Real & _F_RT;
  const Real & _Vm_RT;

};

