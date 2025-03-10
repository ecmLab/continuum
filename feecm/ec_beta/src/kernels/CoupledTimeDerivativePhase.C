#include "CoupledTimeDerivativePhase.h"
#include "Function.h"
registerADMooseObject("ecBetaApp", CoupledTimeDerivativePhase);
InputParameters
CoupledTimeDerivativePhase::validParams()
{
	InputParameters params= ADKernel::validParams();
	params.addClassDescription("Coupled Time Derivative for Poisson Equation");
	params.addParam<Real>("n",1,"Electron Transfer");
	params.addParam<Real>("F",1,"Faraday Constant");
	params.addParam<MaterialPropertyName>("C",1,"Equilibrium Coefficient");
	params.addParam<Real>("scale",1,"Scale Factor");
	params.addRequiredCoupledVar("phi","Order Parameter");
	return params;
}
CoupledTimeDerivativePhase::CoupledTimeDerivativePhase(const InputParameters & parameters) : ADKernel(parameters),
_n(getParam<Real>("n")),
_F(getParam<Real>("F")),
_scale(getParam<Real>("scale")),
_C(getADMaterialProperty<Real>("C")),
_phi_var_dot(adCoupledDot("phi"))
{
}
ADReal
CoupledTimeDerivativePhase::computeQpResidual()
{
	return _n*_F*_scale*_C[_qp]*_phi_var_dot[_qp] * _test[_i][_qp];
}
