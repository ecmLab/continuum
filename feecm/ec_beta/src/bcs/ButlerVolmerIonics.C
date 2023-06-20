
#include "ButlerVolmerIonics.h"

registerADMooseObject("ecBetaApp", ButlerVolmerIonics);

InputParameters
ButlerVolmerIonics::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addClassDescription("Butler-Volmer boundary condition of the Ionics material at the Interface with the reference material.");

// Add a parameter with a default value; this value can be overridden in the input file.
  params.addParam<Real>("F_RT", 0.03868, "The constant of F/RT when T = 300K, in unit 1/mV");
  params.addParam<Real>("reaction_rate", 0.5, "The reaction_rate of Li+ to Li, unitless.");
  params.addParam<Real>("LiPotRef", 0.0, "The Li potental in the reference material, in mV. By default the reference material is Li-Metal");
  params.addRequiredParam<Real>("ex_current", "The exchange current density at the interface, in unit mA/cm^2.");

  return params;
}

ButlerVolmerIonics::ButlerVolmerIonics(const InputParameters & parameters)
  : ADIntegratedBC(parameters),

  // Get the parameters from the input file
   _F_RT(getParam<Real>("F_RT")),
   _reaction_rate(getParam<Real>("reaction_rate")),
   _LiPotRef(getParam<Real>("LiPotRef")),
   _exchange_current(getParam<Real>("ex_current"))

{
}

ADReal
ButlerVolmerIonics::computeQpResidual()
{
  ADReal k1 = std::exp(_reaction_rate * _F_RT * (_LiPotRef - _u[_qp]));
  ADReal k2 = std::exp(- (1 - _reaction_rate) * _F_RT * (_LiPotRef - _u[_qp]));
  return -_test[_i][_qp] * _exchange_current * (k1 - k2);
//  if (k1<k2) {
//    return -_test[_i][_qp] * _exchange_current * (k1 - k2);
//  } else {
//    return -_test[_i][_qp] * 0.0;
//  }
}
