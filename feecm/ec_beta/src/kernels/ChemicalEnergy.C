#include "ChemicalEnergy.h"

registerADMooseObject("ecBetaApp", ChemicalEnergy);
InputParameters
ChemicalEnergy::validParams()
{
	InputParameters params = ADKernel::validParams();
	params.addClassDescription("Computes the Chemical Free Energy");
	params.addRequiredParam<MaterialPropertyName>("w","the coefficient associated with chemical composition");
	params.addParam<Real>("scale",1, "scaling parameter");
	return params;
}
ChemicalEnergy::ChemicalEnergy(const InputParameters & parameters) : ADKernel(parameters),
_w(getADMaterialProperty<Real>("w")),
_scale(getParam<Real>("scale"))
{
}
ADReal
ChemicalEnergy::computeQpResidual()
{
	return _test[_i][_qp]*_scale*_w[_qp]*30*std::pow(-1 + _u[_qp],2)*std::pow(_u[_qp],2);
}
