
#pragma once

#include "ADKernel.h"

//class ADNernstPlanckConvection;

class ADMassFluxEEL : public ADKernel
{
public:
  static InputParameters validParams();
  
  ADMassFluxEEL(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual();
  
  const ADVariableGradient & _grad_V;

  const MaterialProperty<Real> & _conductivity;
  const MaterialProperty<Real> & _mobility;
  const Real & _F;
  const Real & _R;
  const Real & _T;
  const Real & _c0;
  const Real & _scale;
};

