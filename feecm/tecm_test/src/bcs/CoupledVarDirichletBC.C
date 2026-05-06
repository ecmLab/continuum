#include "CoupledVarDirichletBC.h"

registerMooseObject("tecm_testApp", CoupledVarDirichletBC);

InputParameters
CoupledVarDirichletBC::validParams()
{
  InputParameters params = DirichletBCBase::validParams();
  params.addRequiredCoupledVar("value", "Value of the BC");
  params.addClassDescription("Imposes the essential boundary condition $u=g$");
  return params;
}

CoupledVarDirichletBC::CoupledVarDirichletBC(const InputParameters & parameters)
  : DirichletBCBase(parameters), _value(coupledValue("value"))
{
}

Real
CoupledVarDirichletBC::computeQpValue()
{
  return _value[_qp];
}
