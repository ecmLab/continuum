// Copyright 2025, CEWLAB, All Rights Reserved
// License: L-GPL 3.0
#pragma once

#include "ADKernelValue.h"

class MaterialSource : public ADKernelValue
{
public:
  static InputParameters validParams();

  MaterialSource(const InputParameters & parameters);

protected:
  virtual ADReal precomputeQpResidual() override;

  const ADMaterialProperty<Real> & _prop;

  const Real _coef;
};
