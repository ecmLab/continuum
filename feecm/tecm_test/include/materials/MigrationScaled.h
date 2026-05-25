// Copyright 2025, CEWLAB, All Rights Reserved
#pragma once

#include "ElectroChemicalEnergyDensity.h"

class MigrationScaled : public ElectroChemicalEnergyDensity
{
public:
  static InputParameters validParams();

  MigrationScaled(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  /// The electric conductivity
  const ADMaterialProperty<Real> & _sigma;
  const Real _F;
  const Real _R;
  const Real _T;
  const ADMaterialProperty<Real> & _L0;
  const ADMaterialProperty<Real> & _phi0;
};
