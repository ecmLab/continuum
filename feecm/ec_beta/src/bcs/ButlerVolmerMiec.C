
#include "ButlerVolmerMiec.h"

registerADMooseObject("ecBetaApp", ButlerVolmerMiec);

InputParameters
ButlerVolmerMiec::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addClassDescription("Butler-Volmer boundary condition of the MIEC material at the Interface with the reference material.");

// Add a parameter with a default value; this value can be overridden in the input file.
  params.addParam<Real>("F_RT", 0.03868, "The constant of F/RT when T = 300K, in unit 1/mV");
  params.addParam<Real>("reaction_rate", 0.5, "The reaction_rate of Li+ to Li, unitless.");
  params.addParam<Real>("LiCrtRef", 0.0, "The current density in the reference material, in mA/cm^2. It is possible external current is applied on the reference electrode material. By default no external current exist");
  params.addParam<Real>("LiPotRef", 0.0, "The Li potental in the reference material, in mV. By default the other material is Li-Metal");
  params.addRequiredParam<Real>("ex_current", "The exchange current density at the interface, in unit mA/cm^2.");

// Add a coupled parameter: potEn
  params.addRequiredCoupledVar("potEn", "The potential of Electron");

  return params;
}

ButlerVolmerMiec::ButlerVolmerMiec(const InputParameters & parameters)
  : ADIntegratedBC(parameters),

  // Get the parameters from the input file
   _F_RT(getParam<Real>("F_RT")),
   _reaction_rate(getParam<Real>("reaction_rate")),
   _LiCrtRef(getParam<Real>("LiCrtRef")),
   _LiPotRef(getParam<Real>("LiPotRef")),
   _exchange_current(getParam<Real>("ex_current")),
   
  // Couple to the potential of electron
   _potEn(adCoupledValue("potEn"))
{
}

ADReal
ButlerVolmerMiec::computeQpResidual()
{
  ADReal k1 = std::exp(_reaction_rate * _F_RT * (_LiPotRef - _u[_qp] - _potEn[_qp]));
  ADReal k2 = std::exp(- (1 - _reaction_rate) * _F_RT * (_LiPotRef - _u[_qp] - _potEn[_qp]));
  if (k1<k2) {
    return _test[_i][_qp] * (_LiCrtRef - _exchange_current * (k1 - k2));
  } else {
    return _test[_i][_qp] * _LiCrtRef;
  }
}
