// Copyright 2025, CEWLAB, All Rights Reserved
#include "LithiationDeformationGradient.h"

registerADMooseObject("tecm_testApp", LithiationDeformationGradient);

InputParameters
LithiationDeformationGradient::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Computes the eigen deformation gradient associated with isotropic "
      "lithiation-induced volumetric breathing: "
      "Jl = 1 + alpha_Li * soc, Fl = cbrt(Jl) * I. "
      "Sign convention: alpha_Li > 0 for materials that expand on lithiation "
      "(e.g. graphite anode), alpha_Li < 0 for materials that contract on "
      "lithiation (or expand on delithiation).");
  params.addRequiredParam<MaterialPropertyName>("lithiation_deformation_gradient",
                                                "Name of the lithiation deformation gradient");
  params.addRequiredCoupledVar("state_of_charge",
                               "Fractional state of charge in [0, 1]");
  params.addRequiredParam<MaterialPropertyName>(
      "lithiation_expansion_coefficient",
      "Volumetric expansion coefficient alpha_Li (block-dependent material property)");
  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

LithiationDeformationGradient::LithiationDeformationGradient(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _Fl_name(getParam<MaterialPropertyName>("lithiation_deformation_gradient")),
    _Fl(declareADPropertyByName<RankTwoTensor>(_Fl_name)),
    _soc(adCoupledValue("state_of_charge")),
    _alpha_Li(getADMaterialProperty<Real>("lithiation_expansion_coefficient"))
{
}

void
LithiationDeformationGradient::computeQpProperties()
{
  const ADReal Jl = 1.0 + _alpha_Li[_qp] * _soc[_qp];
  _Fl[_qp] = std::cbrt(Jl) * ADRankTwoTensor::Identity();
}
