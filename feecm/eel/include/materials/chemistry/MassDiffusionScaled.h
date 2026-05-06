// Copyright 2023, UChicago Argonne, LLC All Rights Reserved
// License: L-GPL 3.0
#pragma once

#include "Material.h"
#include "DualChemicalEnergyDensity.h"

class MassDiffusionScaled : public DualChemicalEnergyDensity
{
public:
  static InputParameters validParams();

  MassDiffusionScaled(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  /// The mobility
  const ADMaterialProperty<Real> & _M;
  const ADMaterialProperty<Real> &  _phi0;
  const ADMaterialProperty<Real> &  _L0;
};
