#include "FBulk.h"
#include "libmesh/utility.h"
#include <vector>

registerMooseObject("ecBetaApp", FBulk);

InputParameters FBulk::validParams()
{

  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates an integral over the bulk");
  params.addRequiredCoupledVar("op", "The order Parameter");
  params.addParam<MaterialPropertyName>("A", "The coefficients of the energy");
  params.addParam<Real>("len_scale",1.0,"the len_scale of the unit");
  return params;
}

FBulk::FBulk(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  _op(coupledValue("op")),
  _A(getMaterialProperty<Real>("A")),
  _len_scale(getParam<Real>("len_scale"))
{
}

Real
FBulk::computeQpIntegral()
{
  return -_A[_qp]*std::pow(_op[_qp],2)*std::pow((1 - _op[_qp]),2);
}
