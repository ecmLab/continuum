// Copyright 2025, ToBeDecided, All Rights Reserved
// License: L-GPL 3.0
#pragma once

#include "Material.h"
#include "DerivativeMaterialInterface.h"
#include "TecmUtils.h"
#include <Eigen/Dense>

class ChemicalPotential : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  ChemicalPotential(const InputParameters & parameters);

  void computeProperties() override;

protected:
  /// Chemical potential
  ADMaterialProperty<Real> & _mu;

  /// Chemical potential gradient
  ADMaterialProperty<RealVectorValue> & _grad_mu;

  /// Energy names
  const std::vector<MaterialPropertyName> _psi_names;

  /// Energy derivatives
  std::vector<const ADMaterialProperty<Real> *> _d_psi_d_c_dot;

  /// Concentration
  const MooseVariable * _c_var;

  /// the current test function
  const VariableTestValue & _test;

  /// gradient of the test function
  const VariableTestGradient & _grad_test;

private:
  TecmUtils::ADRealEigenVector L2Projection();
};
