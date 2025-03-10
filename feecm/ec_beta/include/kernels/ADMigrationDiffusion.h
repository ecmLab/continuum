
// Including the "Diffusion" Kernel here so we can extend it
#include "ADKernel.h"

class ADMigrationDiffusion : public ADKernel
{
public:
  static InputParameters validParams();

  ADMigrationDiffusion(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  //MaterialPropertyName _diffusivity;
  //const ADMaterialProperty<Real> & _diffusivity_coef;
  const ADVariableValue & _c;
  const ADVariableGradient & _grad_c;
  //const ADVariableValue & _T;

  const MaterialProperty<Real> & _conductivity;
  const Real & _c0;
  const Real & _zIons;
  //const Real & _R;
  const Real & _F_RT;
  const Real & _scale;
};
