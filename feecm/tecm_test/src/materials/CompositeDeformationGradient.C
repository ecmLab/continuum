// Copyright 2025, CEWLAB, All Rights Reserved
#include "CompositeDeformationGradient.h"

registerADMooseObject("tecm_testApp", CompositeDeformationGradient);

InputParameters
CompositeDeformationGradient::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Computes a composite deformation gradient as the ordered product of a "
      "list of input deformation gradients. Useful for composing multiple "
      "eigen deformation gradients (e.g. cure shrink * Li breathing) into a "
      "single property to feed into MechanicalDeformationGradient.");
  params.addRequiredParam<MaterialPropertyName>(
      "composite_deformation_gradient",
      "Name of the output composite deformation gradient");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "factors",
      "Names of the input deformation gradient material properties to multiply, in order");
  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

CompositeDeformationGradient::CompositeDeformationGradient(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _F_name(getParam<MaterialPropertyName>("composite_deformation_gradient")),
    _F(declareADPropertyByName<RankTwoTensor>(_F_name))
{
  for (const auto & nm : getParam<std::vector<MaterialPropertyName>>("factors"))
    _factors.push_back(&getADMaterialProperty<RankTwoTensor>(nm));
}

void
CompositeDeformationGradient::computeQpProperties()
{
  _F[_qp].setToIdentity();
  for (const auto * Fi : _factors)
    _F[_qp] *= (*Fi)[_qp];
}
