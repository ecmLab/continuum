// Created by: Zeeshan Ahmad

#include "ChemPotValue.h"

registerMooseObject("edl_solidApp", ChemPotValue);

InputParameters
ChemPotValue::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Kernel to set the value of mu pointwise equal to a target value.");
  params.addRequiredParam<MaterialPropertyName>("mutarget",
                                                "The target value of the chemical potential");
  params.addRequiredParam<MaterialPropertyName>(
      "mu", "The name of the chemical potential to set equal to the target value");
  params.addCoupledVar("args", "Vector of nonlinear variable arguments this object depends on");
  return params;
}

ChemPotValue::ChemPotValue(const InputParameters & parameters)
  : ADKernel(parameters),
    _mutarget(getADMaterialProperty<Real>("mutarget")),
    _mu(getADMaterialProperty<Real>("mu"))
{
}

ADReal
ChemPotValue::computeQpResidual()
{
  return (_mu[_qp] - _mutarget[_qp]) * _test[_i][_qp];
}