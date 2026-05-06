// Created by: Zeeshan Ahmad

#pragma once

#include "ADKernel.h"
#include "DerivativeMaterialPropertyNameInterface.h"

class ChemPotValueHC : public ADKernel, public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  ChemPotValueHC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual();

  /// Mobility
  const Real & _mutarget;
  const ADMaterialProperty<Real> & _mu0;
  const ADVariableValue &_cb, &_phi;
  const ADMaterialProperty<Real> & _z;
  const Real & _kT;

  //   /// Interfacial parameter
  //   const ADMaterialProperty<Real> & _kappa;
};