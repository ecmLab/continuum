#include "NucleationEnergy.h"
registerADMooseObject("ecBetaApp", NucleationEnergy);
InputParameters
NucleationEnergy::validParams()
{
	InputParameters params=ADKernel::validParams();
	params.addClassDescription("Computing Nucleation Energy");
	params.addRequiredParam<MaterialPropertyName>("delta_G", "Nucleation Energy Coefficient");
	params.addRequiredParam<MaterialPropertyName>("gamma","Penalty Factor for Nucleation");
	params.addParam<Real>("scale",1,"Scale factor");
	return params;
}
NucleationEnergy::NucleationEnergy(const InputParameters & parameters) : ADKernel(parameters),
_delta_G(getADMaterialProperty<Real>("delta_G")),
_gamma(getADMaterialProperty<Real>("gamma")),
_scale(getParam<Real>("scale"))
{
}
ADReal
NucleationEnergy::computeQpResidual()
{
	return 2*_test[_i][_qp]*_scale*_delta_G[_qp]*_gamma[_qp]*(_u[_qp]-1);
}
