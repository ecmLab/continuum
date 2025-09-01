// Copyright 2025, ToBeDecided, All Rights Reserved
// License: L-GPL 3.0
#pragma once

#include "ThermodynamicForce.h"

class FirstPiolaKirchhoffStress : public ThermodynamicForce<RankTwoTensor>
{
public:
  static InputParameters validParams();

  FirstPiolaKirchhoffStress(const InputParameters & parameters);
};
