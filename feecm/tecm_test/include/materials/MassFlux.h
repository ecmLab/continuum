// Copyright 2025, CEWLAB, All Rights Reserved
// License: L-GPL 3.0
#pragma once

#include "ThermodynamicForce.h"

class MassFlux : public ThermodynamicForce<RealVectorValue>
{
public:
  static InputParameters validParams();

  MassFlux(const InputParameters & parameters);
};
