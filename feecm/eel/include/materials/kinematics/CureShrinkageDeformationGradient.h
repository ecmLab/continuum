// Copyright 2023, UChicago Argonne, LLC All Rights Reserved
// License: L-GPL 3.0
#pragma once

#include "Material.h"
#include "ADRankTwoTensorForward.h"
#include "DerivativeMaterialInterface.h"

class CureShrinkageDeformationGradient : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  CureShrinkageDeformationGradient(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Name of the cure-shrinkage deformation gradient
  const MaterialPropertyName _Fc_name;

  /// The cure-shrinkage deformation gradient
  ADMaterialProperty<RankTwoTensor> & _Fc;

  /// Cure variable, expected to range in [0, 1]
  const ADVariableValue & _c;

  /// Volumetric shrinkage at full cure (positive; e.g. 0.05 for 5%), block-dependent
  const ADMaterialProperty<Real> & _eps_sh;
};
