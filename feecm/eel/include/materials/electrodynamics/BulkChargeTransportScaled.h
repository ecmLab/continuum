// Copyright 2023, UChicago Argonne, LLC All Rights Reserved
// License: L-GPL 3.0
#pragma once

#include "ElectricalEnergyDensity.h"

class BulkChargeTransportScaled : public ElectricalEnergyDensity
{
public:
  static InputParameters validParams();

  BulkChargeTransportScaled(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  /// The electric conductivity
  const ADMaterialProperty<Real> & _sigma;
  const ADMaterialProperty<Real> & _phi0;
  const ADMaterialProperty<Real> & _L0;

  /// The derivative of the electrical energy density w.r.t. the log temperature
  ADMaterialProperty<Real> * _d_E_d_lnT;

};
