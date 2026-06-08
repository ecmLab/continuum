// Copyright 2025, CEWLAB, All Rights Reserved
#pragma once

#include "ThermodynamicForce.h"

class MassSource : public ThermodynamicForce<Real>
{
public:
  static InputParameters validParams();

  MassSource(const InputParameters & parameters);
};
