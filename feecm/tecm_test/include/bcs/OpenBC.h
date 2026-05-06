// Copyright 2025, CEWLAB, All Rights Reserved
// License: L-GPL 3.0
#pragma once

#include "ADIntegratedBC.h"

class OpenBC : public ADIntegratedBC
{
public:
  OpenBC(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

  const ADMaterialProperty<RealVectorValue> & _flux;
};
