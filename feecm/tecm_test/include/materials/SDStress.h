// Copyright 2025, CEWLAB, All Rights Reserved
#pragma once

#include "ThermodynamicForce.h"

class SDStress : public ThermodynamicForce<RankTwoTensor>
{
public:
  static InputParameters validParams();

  SDStress(const InputParameters & parameters);
};
