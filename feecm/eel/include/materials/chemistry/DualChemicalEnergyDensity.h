// Copyright 2023, UChicago Argonne, LLC All Rights Reserved
// License: L-GPL 3.0
#pragma once

#include "Material.h"
#include "DerivativeMaterialInterface.h"

class DualChemicalEnergyDensity : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  DualChemicalEnergyDensity(const InputParameters & parameters);

protected:
  /// Name of the dual chemical energy density
  const MaterialPropertyName _energy_name;

  /// The chemical potential
  const MaterialPropertyName _mu_name;

  /// The gradient of the chemical potential
  const ADMaterialProperty<RealVectorValue> & _grad_mu;

  /// The dual chemical energy density
  ADMaterialProperty<Real> & _zeta;

  /// Derivative of the dual chemical energy density w.r.t. the chemical potential gradient
  ADMaterialProperty<RealVectorValue> & _d_zeta_d_grad_mu;
};
