// Created by: Zeeshan Ahmad

#pragma once

#include "ADKernel.h"

class ChemPotValue : public ADKernel, public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  ChemPotValue(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual();

  const ADMaterialProperty<Real> & _mutarget;
  const ADMaterialProperty<Real> & _mu;
};