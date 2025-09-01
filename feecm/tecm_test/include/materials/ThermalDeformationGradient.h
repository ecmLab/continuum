// Copyright 2023, UChicago Argonne, LLC All Rights Reserved
// License: L-GPL 3.0
#pragma once

#include "Material.h"
#include "ADRankTwoTensorForward.h"
#include "Function.h"
#include "DerivativeMaterialInterface.h"

class ThermalDeformationGradient : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  ThermalDeformationGradient(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Name of the thermal deformation gradient
  const MaterialPropertyName _Ft_name;

  // The thermal deformation gradient
  ADMaterialProperty<RankTwoTensor> & _Ft;

  /// Temperature variable name
  const VariableName _T_name;

  // The current temperature
  const ADVariableValue & _T;

  // The reference temperature
  const VariableValue & _T_ref;

  // The thermal expansion coefficient
  const ADMaterialProperty<Real> & _alpha_t;

  /// Derivative of the thermal deformation gradient jacobian w.r.t. the temperature
  ADMaterialProperty<Real> & _d_Jt_d_T;
};
