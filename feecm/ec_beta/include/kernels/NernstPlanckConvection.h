
#pragma once

#include "Kernel.h"

class NernstPlanckConvection;

class NernstPlanckConvection : public Kernel
{
public:
  static InputParameters validParams();
  
  NernstPlanckConvection(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// int label for temperature variable
  unsigned int _V_var;

  /// Coupled variable for the temperature
  const VariableValue & _V;

  /// Variable gradient for temperature
  const VariableGradient & _grad_V;

  /// Diffusivity material property
  const MaterialProperty<Real> & _diffusivity;

  /// Will be set from the input file
  const Real & _zIons;
  const Real & _F_RT;

  //Real _Qh;
  //Real _kb;
};

