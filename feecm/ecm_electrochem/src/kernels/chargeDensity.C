#include "chargeDensity.h"

registerADMooseObject("liExpulsionApp", chargeDensity);

InputParameters
chargeDensity::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Implements the weak form $(\\psi_i, -f)$.");

// Add a coupled parameter: Concentration of the ions
  params.addRequiredCoupledVar("conIons", "The concentration of the ions");

    params.addParam<Real>("zIons", 1, "The charge state of the ions, default be positive 1");
    params.addParam<Real>("scale", 1, "The scale of the equation, default no scaling");

  return params;
}

chargeDensity::chargeDensity(const InputParameters & parameters)
  : ADKernel(parameters),
   _conIons(adCoupledValue("conIons")),
   // Couple to the gradient of the concentration
  // _grad_con(coupledGradient("conIons")),
   _zIons(getParam<Real>("zIons")),
   _scale(getParam<Real>("scale"))
{
}

ADReal
chargeDensity::computeQpResidual()
{
 return -_scale * _test[_i][_qp] * _zIons * _conIons[_qp];
}

//Real
//chargeDensity::computeQpJacobian()
//{
//  return -_scale * _test[_i][_qp] * _zIons * _grad_con[_qp]);
//}
