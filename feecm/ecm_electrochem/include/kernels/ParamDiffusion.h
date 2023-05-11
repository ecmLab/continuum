
#pragma once

// Including the "Diffusion" Kernel here so we can extend it
#include "ADKernel.h"
class paramDiffusion;

class paramDiffusion : public ADKernel
{
public:
  static InputParameters validParams();

  paramDiffusion(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  MaterialPropertyName _diffusivity;
  const ADMaterialProperty<Real> & _diffusivity_coef;

//  usingKernelMembers;
};

