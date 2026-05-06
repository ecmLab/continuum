// Copyright 2023, UChicago Argonne, LLC All Rights Reserved
// License: L-GPL 3.0
#pragma once

#include "Material.h"
#include "ADRankTwoTensorForward.h"
#include "DerivativeMaterialInterface.h"

class LithiationDeformationGradient : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  LithiationDeformationGradient(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Name of the lithiation deformation gradient
  const MaterialPropertyName _Fl_name;

  /// The lithiation deformation gradient
  ADMaterialProperty<RankTwoTensor> & _Fl;

  /// State of charge in [0, 1]
  const ADVariableValue & _soc;

  /// Volumetric expansion coefficient on lithiation (positive expands, negative shrinks),
  /// block-dependent
  const ADMaterialProperty<Real> & _alpha_Li;
};
