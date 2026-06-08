#pragma once

#include "ElementIntegralPostprocessor.h"

class FBulk: public ElementIntegralPostprocessor
{
public:
  FBulk(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpIntegral();
  const VariableValue & _op;
  const MaterialProperty<Real> & _A;
  const Real _len_scale;

};