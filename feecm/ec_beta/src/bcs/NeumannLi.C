
#include "NeumannLi.h"

registerADMooseObject("ecBetaApp", NeumannLi);

InputParameters
NeumannLi::validParams()
{

  InputParameters params = ADIntegratedBC::validParams();
  params.addClassDescription("The Neumann boundary condition for Li+ at Cathode side.");

  params.addRequiredCoupledVar("potEn","The variable representing the potential of electron.");

  return params;
}

NeumannLi::NeumannLi(const InputParameters & parameters)
  : ADIntegratedBC(parameters),

// Get the gradient of the variable
    _potEn_grad(adCoupledGradient("potEn")),

// Get the parameters from the input file
  _electronic_conductivity(getADMaterialProperty<Real>("electronic_conductivity")),
  _inlet_current(getADMaterialProperty<Real>("inlet_current"))
{
}

ADReal
NeumannLi::computeQpResidual()
{
  return _test[_i][_qp] * (10000 * _electronic_conductivity[_qp] * _grad_u[_qp] * _normals[_qp] - _inlet_current[_qp]);
}
