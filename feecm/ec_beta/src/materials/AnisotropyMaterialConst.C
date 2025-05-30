#include "AnisotropyMaterialConst.h"

registerMooseObject("ecBetaApp", AnisotropyMaterialConst);


InputParameters 
AnisotropyMaterialConst::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Anisotropy Material Const");
  params.addRequiredParam<Real>("sig11", "Value for the (1,1) component of the tensor");
  params.addRequiredParam<Real>("sig12", "Value for the (1,2) component of the tensor");
  params.addRequiredParam<Real>("sig13", "Value for the (1,3) component of the tensor");
  params.addRequiredParam<Real>("sig21", "Value for the (2,1) component of the tensor");
  params.addRequiredParam<Real>("sig22", "Value for the (2,2) component of the tensor");
  params.addRequiredParam<Real>("sig23", "Value for the (2,3) component of the tensor");
  params.addRequiredParam<Real>("sig31", "Value for the (3,1) component of the tensor");
  params.addRequiredParam<Real>("sig32", "Value for the (3,2) component of the tensor");
  params.addRequiredParam<Real>("sig33", "Value for the (3,3) component of the tensor");

  return params;
}
AnisotropyMaterialConst::AnisotropyMaterialConst(const InputParameters & parameters)
  : Material(parameters),
    _sig11(getParam<Real>("sig11")),
    _sig12(getParam<Real>("sig12")),
    _sig13(getParam<Real>("sig13")),
    _sig21(getParam<Real>("sig21")),
    _sig22(getParam<Real>("sig22")),
    _sig23(getParam<Real>("sig23")),
    _sig31(getParam<Real>("sig31")),
    _sig32(getParam<Real>("sig32")),
    _sig33(getParam<Real>("sig33")),
    _tensor_property(declareADProperty<RealTensorValue>("tensor_property"))

{
}

void AnisotropyMaterialConst::computeQpProperties()
{
  // Initialize a RankTwoTensor with given values
  RankTwoTensor local_tensor;
  local_tensor(0, 0) = _sig11; // (1,1) component
  local_tensor(0, 1) = _sig12; // (1,2) component
  local_tensor(0, 2) = _sig13; // (1,3) component
  local_tensor(1, 0) = _sig21; // (2,1) component
  local_tensor(1, 1) = _sig22; // (2,2) component
  local_tensor(1, 2) = _sig23; // (2,3) component
  local_tensor(2, 0) = _sig31; // (3,1) component
  local_tensor(2, 1) = _sig32; // (3,2) component
  local_tensor(2, 2) = _sig33; // (3,3) component
  

  // Assign to material property
  _tensor_property[_qp] = local_tensor;
}
