
#include "ButlerVolmerMiecSrf.h"

registerADMooseObject("ecBetaApp", ButlerVolmerMiecSrf);

InputParameters
ButlerVolmerMiecSrf::validParams()
{
  InputParameters params = ADIntegratedBC::validParams(); 
  params.addClassDescription("Butler-Volmer boundary condition of the MIEC material at the surface. A better version need to be derived in the future");

// Add a parameter with a default value; this value can be overridden in the input file.
  params.addParam<Real>("F_RT", 0.03868, "The constant of F/RT when T = 300K, in unit 1/mV");
  params.addParam<Real>("reaction_rate", 0.5, "The reaction_rate of Li+ to Li, unitless.");
  params.addParam<Real>("electron_concentration", 0.008, "The electron concentration inside the miec, unitless.");
  params.addRequiredParam<Real>("ex_current", "The exchange current density at the interface, in unit mA/cm^2.");

// Add a coupled parameter: potEn
  params.addRequiredCoupledVar("potEn", "The potential of Electron");

  return params;
}

ButlerVolmerMiecSrf::ButlerVolmerMiecSrf(const InputParameters & parameters)
  : ADIntegratedBC(parameters),

  // Get the parameters from the input file
   _F_RT(getParam<Real>("F_RT")),
   _reaction_rate(getParam<Real>("reaction_rate")),
   _ele_conc(getParam<Real>("electron_concentration")),
   _exchange_current(getParam<Real>("ex_current")),

  // Couple to the potential of electron
   _potEn(adCoupledValue("potEn"))
{
}

ADReal
ButlerVolmerMiecSrf::computeQpResidual()
{
  ADReal k1 = std::exp(_reaction_rate * _F_RT * (-_u[_qp] - _potEn[_qp]));
  ADReal k2 = std::exp(- (1 - _reaction_rate) * _F_RT * (-_u[_qp] - _potEn[_qp]));

  if (k1<k2) {
    return -_test[_i][_qp] * _ele_conc * _exchange_current * (k1-k2);
  } else {
    return -_test[_i][_qp] * 0.0;
  }
}
