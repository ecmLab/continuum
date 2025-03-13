// Copyright 2023, UChicago Argonne, LLC All Rights Reserved
// License: L-GPL 3.0
#pragma once

#include "MechanicalEnergyDensity.h"

class NeoHookeanSolid : public MechanicalEnergyDensity
{
public:
  static InputParameters validParams();

  NeoHookeanSolid(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  /// Lame's first parameter
  const ADMaterialProperty<Real> & _lambda;

  /// Shear modulus
  const ADMaterialProperty<Real> & _G;

  /// Non-swelling portion of the pressure
  ADMaterialProperty<Real> * _p;
};
