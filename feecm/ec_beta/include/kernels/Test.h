// Including the "Diffusion" Kernel here so we can extend it
#include "ADKernel.h"

class Test : public ADKernel
{
public:
  static InputParameters validParams();

  Test(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const Real & _A;

};
