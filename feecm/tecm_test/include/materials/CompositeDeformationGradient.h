// Copyright 2025, CEWLAB, All Rights Reserved
#pragma once

#include "Material.h"
#include "ADRankTwoTensorForward.h"
#include "DerivativeMaterialInterface.h"

class CompositeDeformationGradient : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  CompositeDeformationGradient(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Name of the output composite deformation gradient
  const MaterialPropertyName _F_name;

  /// The composite deformation gradient
  ADMaterialProperty<RankTwoTensor> & _F;

  /// Input deformation gradients to multiply, in order
  std::vector<const ADMaterialProperty<RankTwoTensor> *> _factors;
};
