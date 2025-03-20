
// Including the "Diffusion" Kernel here so we can extend it
#include "ADKernel.h"

class ADCurrentFluxEEL : public ADKernel
{
public:
  static InputParameters validParams();

  ADCurrentFluxEEL(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const ADVariableValue & _c;
  const ADVariableGradient & _grad_c;
  const MaterialProperty<Real> & _conductivity;
  const Real & _F;
  const Real & _R;
  const Real & _T;
  const Real & _c0;
  const Real & _scale;
};
