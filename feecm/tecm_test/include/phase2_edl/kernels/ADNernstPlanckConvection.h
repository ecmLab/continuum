#pragma once

#include "ADKernel.h"

/**
 * Automatic-differentiation kernel implementing the electromigration term
 * ∇·( (z_i F D_i)/(R T) c_i ∇φ ) for the TECM EDL model.
 */
class ADNernstPlanckConvection : public ADKernel
{
public:
  static InputParameters validParams();

  ADNernstPlanckConvection(const InputParameters & parameters);

protected:
  ADReal computeQpResidual() override;

private:
  const ADMaterialProperty<Real> & _diffusivity;
  const ADVariableGradient & _grad_potential;
  const ADVariableValue & _temperature;
  const Real _valence;
  const Real _faraday_constant;
  const Real _gas_constant;
};
