#include "ElectrodeDrivingForce.h"

registerADMooseObject("ecBetaApp", ElectrodeDrivingForce);
InputParameters
ElectrodeDrivingForce::validParams()
{
	InputParameters params = ADKernel::validParams();
	params.addClassDescription("Computes the Electrode Driving Force Energy");
	params.addRequiredParam<MaterialPropertyName>("h","the coefficient associated with chemical composition");
	params.addParam<Real>("alpha",1, "reaction coefficient");
	params.addParam<Real>("beta",1, "reaction coefficient");
	params.addParam<Real>("n",1, "Charge Transfer Number");
	params.addParam<Real>("R",8.84, "Gas Constant");
	params.addParam<Real>("T",298, "Temperature");
	params.addParam<Real>("F",96485, "Faraday Constant");
	params.addRequiredCoupledVar("pot","Potential Field");
	params.addRequiredCoupledVar("ref_pot","Ref Potential Field");
	params.addRequiredCoupledVar("conc","Concentration Field");
	params.addParam<Real>("scale",1, "scaling parameter");
	return params;
}
ElectrodeDrivingForce::ElectrodeDrivingForce(const InputParameters & parameters) : ADKernel(parameters),
_h(getADMaterialProperty<Real>("h")),
_alpha(getParam<Real>("alpha")),
_beta(getParam<Real>("beta")),
_n(getParam<Real>("n")),
_R(getParam<Real>("R")),
_T(getParam<Real>("T")),
_F(getParam<Real>("F")),
_pot(coupledValue("pot")),
_ref_pot(coupledValue("ref_pot")),
_conc(coupledValue("conc")),
_scale(getParam<Real>("scale"))
{
}
ADReal
ElectrodeDrivingForce::computeQpResidual()
{
	Real dG = _n*_F*(_pot[_qp]-_ref_pot[_qp])/(_R*_T);
	return _test[_i][_qp]*_scale*_h[_qp]*(std::exp((_alpha)*dG)-std::exp(-(_beta)*dG));
}
