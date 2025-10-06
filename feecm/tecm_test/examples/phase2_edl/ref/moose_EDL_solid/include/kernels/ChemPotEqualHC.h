// Created by: Zeeshan Ahmad

#pragma once

#include "ADKernel.h"
#include "DerivativeMaterialPropertyNameInterface.h"

class ChemPotEqualHC : public ADKernel, public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  ChemPotEqualHC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual();

  /// Mobility
  const ADMaterialProperty<Real> &_mu0a, &_mu0b;
  const ADVariableValue &_cb, &_phi;
  const ADMaterialProperty<Real> &_za, &_zb;
  const Real & _kT;

  //   /// Interfacial parameter
  //   const ADMaterialProperty<Real> & _kappa;
};