#pragma once

#include "ADMaterial.h"

// EDLFluxes: provide flux material properties for Poisson and Nernst–Planck
// - elec_flux = epsilon * grad(phi)
// - j_diff_+  = - D_plus  * grad(c_plus)
// - j_diff_-  = - D_minus * grad(c_minus)
// - j_em_+    = (z_plus  F D_plus)/(R T)  * c_plus  * grad(phi)
// - j_em_-    = (z_minus F D_minus)/(R T) * c_minus * grad(phi)
// - charge_density = F (z_plus c_plus + z_minus c_minus)

class EDLFluxes : public ADMaterial
{
public:
  static InputParameters validParams();
  EDLFluxes(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

private:
  // Coupled variables
  const ADVariableGradient & _grad_phi;
  const ADVariableGradient & _grad_c_plus;
  const ADVariableGradient & _grad_c_minus;
  const ADVariableValue & _c_plus_val;
  const ADVariableValue & _c_minus_val;

  // Material inputs
  const ADMaterialProperty<Real> & _permittivity;
  const ADMaterialProperty<Real> & _D_plus;
  const ADMaterialProperty<Real> & _D_minus;

  // Constants
  const Real _F;
  const Real _R;
  const Real _T;
  const Real _z_plus;
  const Real _z_minus;

  // Outputs
  ADMaterialProperty<ADRealVectorValue> & _elec_flux;
  ADMaterialProperty<ADRealVectorValue> & _j_diff_plus;
  ADMaterialProperty<ADRealVectorValue> & _j_diff_minus;
  ADMaterialProperty<ADRealVectorValue> & _j_em_plus;
  ADMaterialProperty<ADRealVectorValue> & _j_em_minus;
  ADMaterialProperty<ADReal> & _charge_density;
};

