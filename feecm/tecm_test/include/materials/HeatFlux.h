// Copyright 2025, CEWLAB, All Rights Reserved
// License: L-GPL 3.0
#pragma once

#include "ThermodynamicForce.h"

class HeatFlux : public ThermodynamicForce<RealVectorValue>
{
public:
  static InputParameters validParams();

  HeatFlux(const InputParameters & parameters);
};
