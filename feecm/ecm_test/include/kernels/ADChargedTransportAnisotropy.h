
// Including the "Diffusion" Kernel here so we can extend it
#include "ADKernel.h"

class ADChargedTransportAnisotropy : public ADKernel
{
public:
  static InputParameters validParams();

  ADChargedTransportAnisotropy(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  MaterialPropertyName _diffusivity;
  const ADMaterialProperty<RealTensorValue> & _diffusivity_coef;

};
