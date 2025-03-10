#include "CouplingFreeEnergy.h"

registerADMooseObject("ecBetaApp", CouplingFreeEnergy);
InputParameters
CouplingFreeEnergy::validParams()
{
	InputParameters params=ADKernel::validParams();
	params.addClassDescription("Evaluates the coupling Free energy");
	params.addRequiredCoupledVar("coupledVar", "The variable from which to compute the gradient component");
	params.addParam<MaterialPropertyName>("c",1,"Coupling Constant");
	params.addParam<Real>("scale",1,"Scaling Parameters");
	return params;
}
CouplingFreeEnergy::CouplingFreeEnergy(const InputParameters & parameters) : ADKernel(parameters),
_coupledVar(coupledValue("coupledVar")),
_c(getADMaterialProperty<Real>("c")),
_scale(getParam<Real>("scale"))
{
}
ADReal
CouplingFreeEnergy::computeQpResidual()
{
	return 2*_scale*_test[_i][_qp]*_c[_qp]*_u[_qp]*std::pow(_coupledVar[_qp],2);
}
