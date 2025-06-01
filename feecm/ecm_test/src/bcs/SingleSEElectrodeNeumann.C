
#include "SingleSEElectrodeNeumann.h"

registerADMooseObject("ecmApp", SingleSEElectrodeNeumann);

InputParameters
SingleSEElectrodeNeumann::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addClassDescription("The Neumann boundary condition for Li+ .");

  return params;
}

SingleSEElectrodeNeumann::SingleSEElectrodeNeumann(const InputParameters & parameters)
  : ADIntegratedBC(parameters),

// Get the parameters from the input file
  _applied_current(getADMaterialProperty<Real>("applied_current"))
{
}

ADReal
SingleSEElectrodeNeumann::computeQpResidual()
{
  return -_test[_i][_qp] * _applied_current[_qp];
}
