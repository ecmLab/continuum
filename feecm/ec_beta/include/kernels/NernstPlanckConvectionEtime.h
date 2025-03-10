
#pragma once

#include "Kernel.h"

class NernstPlanckConvectionEtime;

class NernstPlanckConvectionEtime : public Kernel
{
public:
  static InputParameters validParams();

  NernstPlanckConvectionEtime(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// Coupled variable for the temperature
  const VariableValue & _Efield;

  /// Variable gradient for temperature

  /// Diffusivity material property
  const MaterialProperty<Real> & _diffusivity;

  /// Will be set from the input file
  const Real & _zIons;
  const Real & _F_RT;
  const Real & _scale;

  //Real _Qh;
  //Real _kb;
};
