// Copyright 2023, UChicago Argonne, LLC All Rights Reserved
// License: L-GPL 3.0
#include "CureShrinkageDeformationGradient.h"

registerADMooseObject("EelApp", CureShrinkageDeformationGradient);

InputParameters
CureShrinkageDeformationGradient::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Computes the eigen deformation gradient associated with isotropic "
      "polymerization-induced volumetric shrinkage, parameterized by a cure "
      "variable c in [0,1] and a volumetric shrinkage eps_sh (>=0). "
      "Jc = 1 - c * eps_sh, Fc = cbrt(Jc) * I.");
  params.addRequiredParam<MaterialPropertyName>("cure_shrinkage_deformation_gradient",
                                                "Name of the cure-shrinkage deformation gradient");
  params.addRequiredCoupledVar("cure", "The cure variable, expected to range in [0, 1]");
  params.addRequiredParam<MaterialPropertyName>(
      "volumetric_shrinkage",
      "Material property giving the total volumetric shrinkage at full cure "
      "(e.g. 0.05 for 5%). Allows block-dependent values.");
  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

CureShrinkageDeformationGradient::CureShrinkageDeformationGradient(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _Fc_name(getParam<MaterialPropertyName>("cure_shrinkage_deformation_gradient")),
    _Fc(declareADPropertyByName<RankTwoTensor>(_Fc_name)),
    _c(adCoupledValue("cure")),
    _eps_sh(getADMaterialProperty<Real>("volumetric_shrinkage"))
{
}

void
CureShrinkageDeformationGradient::computeQpProperties()
{
  const ADReal Jc = 1.0 - _c[_qp] * _eps_sh[_qp];
  _Fc[_qp] = std::cbrt(Jc) * ADRankTwoTensor::Identity();
}
