#pragma once

#include "ADInterfaceKernel.h"

// Forward Declarations
class IonicSEInterface : public ADInterfaceKernel
{
public:
  static InputParameters validParams();

  IonicSEInterface(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual(Moose::DGResidualType type) override;
//  virtual ADReal computeQpJacobian(Moose::DGJacobianType type) override;
  
  /// Get parameters from Material system
  const ADMaterialProperty<Real> & _ionic_conductivity;
  const ADMaterialProperty<Real> & _metal_conductivity;
  const ADMaterialProperty<Real> & _exchange_current;
  const ADMaterialProperty<Real> & _reaction_rate;
  const Real & _F_RT;

};

