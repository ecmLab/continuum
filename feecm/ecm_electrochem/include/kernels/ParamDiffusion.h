
#pragma once

// Including the "Diffusion" Kernel here so we can extend it
#include "ADKernel.h"
class ParamDiffusion;

class ParamDiffusion : public ADKernel
{
public:
  static InputParameters validParams();

  ParamDiffusion(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  MaterialPropertyName _conductivity;
  const ADMaterialProperty<Real> & _conductivity_coef;

//  usingKernelMembers;
};

