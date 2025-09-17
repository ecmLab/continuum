#pragma once

#include "ADKernel.h"

/**
 * Nernst-Planck convection kernel for electromigration transport
 * Implements: v · ∇c where v = -μ * z * ∇Ψ = -D * z * F/(RT) * ∇Ψ
 * This is the correct convection form for ion transport in electric fields
 */
class ADNernstPlanckConvection : public ADKernel
{
public:
  static InputParameters validParams();

  ADNernstPlanckConvection(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// Diffusivity material property
  const ADMaterialProperty<Real> & _diffusivity;
  
  /// Ion valence (charge number)
  const Real _valence;
  
  /// Faraday constant
  const Real _faraday_constant;
  
  /// Gas constant
  const Real _gas_constant;
  
  /// Temperature
  const Real _temperature;
  
  /// Electric potential gradient (coupled variable)
  const ADVariableGradient & _potential_gradient;
};