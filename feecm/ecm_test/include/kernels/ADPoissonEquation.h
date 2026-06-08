#include "ADKernel.h"

class ADPoissonEquation : public ADKernel
{
public:
  static InputParameters validParams();
  ADPoissonEquation(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  //const ADVariableGradient & _grad_potential;
  const ADVariableValue & _coupled_concentration;
  const Real _F;
  const Real _sigma;
  const Real _z;
  const Real _scale;
};

