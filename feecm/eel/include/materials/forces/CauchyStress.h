// Copyright 2023, UChicago Argonne, LLC All Rights Reserved
// License: L-GPL 3.0
#pragma once

#include "Material.h"
#include "ADRankTwoTensorForward.h"

class CauchyStress : public Material
{
public:
  static InputParameters validParams();

  CauchyStress(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  ADMaterialProperty<RankTwoTensor> & _cauchy;

  const ADMaterialProperty<RankTwoTensor> & _pk1;

  const ADMaterialProperty<RankTwoTensor> & _F;
};
