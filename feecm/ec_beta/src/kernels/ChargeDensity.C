#include "ChargeDensity.h"

registerADMooseObject("ecBetaApp", ChargeDensity);

InputParameters
ChargeDensity::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Implements the weak form $(\\psi_i, -f)$.");

// Add a coupled parameter: Concentration of the ions
  params.addRequiredCoupledVar("conIons", "The concentration of the ions");

    params.addParam<Real>("zIons", 1, "The Charge state of the ions, default be positive 1");
    params.addParam<Real>("scale", 1, "The scale of the equation, default no scaling");

  return params;
}

ChargeDensity::ChargeDensity(const InputParameters & parameters)
  : ADKernel(parameters),
   _conIons(adCoupledValue("conIons")),
   // Couple to the gradient of the concentration
  // _grad_con(coupledGradient("conIons")),
   _zIons(getParam<Real>("zIons")),
   _scale(getParam<Real>("scale"))
{
}

ADReal
ChargeDensity::computeQpResidual()
{
 return -_scale * _test[_i][_qp] * _zIons * _conIons[_qp];
}

//Real
//ChargeDensity::computeQpJacobian()
//{
//  return -_scale * _test[_i][_qp] * _zIons * _grad_con[_qp]);
//}
