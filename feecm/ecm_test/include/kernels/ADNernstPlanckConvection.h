
#pragma once

#include "ADKernel.h"

//class ADNernstPlanckConvection;

class ADNernstPlanckConvection : public ADKernel
{
public:
  static InputParameters validParams();
  
  ADNernstPlanckConvection(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual();
  
  const ADVariableGradient & _grad_V;
  const MaterialProperty<Real> & _diffusivity;
  const Real & _zIons;
  const Real & _F_RT;
  //const Real & _T;
  const Real & _scale;
};

