// Created by: Zeeshan Ahmad

#include "ChemPotEqual.h"

registerMooseObject("edl_solidApp", ChemPotEqual);

InputParameters
ChemPotEqual::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription(
      "Kernel to enforce pointwise equilibrium of (electro)chemical potentials mua + mub = 0.");
  params.addRequiredParam<MaterialPropertyName>("mua", "The name of the first chemical potential");
  params.addRequiredParam<MaterialPropertyName>("mub", "The name of the second chemical potential");
  params.addCoupledVar("args", "Vector of nonlinear variable arguments this object depends on");
  return params;
}

ChemPotEqual::ChemPotEqual(const InputParameters & parameters)
  : ADKernel(parameters),
    _mua(getADMaterialProperty<Real>("mua")),
    _mub(getADMaterialProperty<Real>("mub"))
{
}

ADReal
ChemPotEqual::computeQpResidual()
{
  return (_mua[_qp] + _mub[_qp]) * _test[_i][_qp];
}
