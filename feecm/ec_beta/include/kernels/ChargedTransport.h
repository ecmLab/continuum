
// Including the "Diffusion" Kernel here so we can extend it
#include "ADKernel.h"

class ChargedTransport : public ADKernel
{
public:
  static InputParameters validParams();

  ChargedTransport(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  MaterialPropertyName _diffusivity;
  const ADMaterialProperty<Real> & _diffusivity_coef;

};
