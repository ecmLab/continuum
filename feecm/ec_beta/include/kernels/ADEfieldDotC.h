#pragma once

#include "ADKernel.h"

class ADEfieldDotC : public ADKernel
{
public:
  static InputParameters validParams();

  ADEfieldDotC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // Coupled variable
  const ADVariableGradient & _grad_phi;

  // Material property
  const ADMaterialProperty<Real> & _coeff;
};
