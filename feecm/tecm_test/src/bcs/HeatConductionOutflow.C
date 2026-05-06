
#include "HeatConductionOutflow.h"

registerMooseObject("tecm_testApp", HeatConductionOutflow);

InputParameters
HeatConductionOutflow::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addClassDescription("Compute the outflow boundary condition.");
  return params;
}

HeatConductionOutflow::HeatConductionOutflow(const InputParameters & parameters)
  : ADIntegratedBC(parameters),
    _thermal_conductivity(getADMaterialProperty<Real>("thermal_conductivity"))
{
}

ADReal
HeatConductionOutflow::computeQpResidual()
{
  return -_test[_i][_qp] * _thermal_conductivity[_qp] * _grad_u[_qp] * _normals[_qp];
}
