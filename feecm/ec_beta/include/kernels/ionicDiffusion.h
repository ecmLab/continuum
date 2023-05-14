
// Including the "Diffusion" Kernel here so we can extend it
#include "ADKernel.h"

class ionicDiffusion : public ADKernel
{
public:
  static InputParameters validParams();

  ionicDiffusion(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  MaterialPropertyName _diffusivity;
  const ADMaterialProperty<Real> & _diffusivity_coef;

};
