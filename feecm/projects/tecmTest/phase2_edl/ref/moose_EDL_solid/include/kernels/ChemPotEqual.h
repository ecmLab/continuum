// Created by: Zeeshan Ahmad

#pragma once

#include "ADKernel.h"
#include "DerivativeMaterialPropertyNameInterface.h"

class ChemPotEqual : public ADKernel, public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  ChemPotEqual(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual();

  const ADMaterialProperty<Real> & _mua;
  const ADMaterialProperty<Real> & _mub;
};