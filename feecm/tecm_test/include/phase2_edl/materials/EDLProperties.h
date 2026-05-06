#pragma once

#include "ADMaterial.h"

/**
 * Provides constant diffusion coefficients and permittivity for the Phase 2 EDL model.
 */
class EDLProperties : public ADMaterial
{
public:
  static InputParameters validParams();
  EDLProperties(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

private:
  const Real _D_plus_value;
  const Real _D_minus_value;
  const Real _permittivity_value;

  ADMaterialProperty<Real> & _D_plus;
  ADMaterialProperty<Real> & _D_minus;
  ADMaterialProperty<Real> & _permittivity;
};

