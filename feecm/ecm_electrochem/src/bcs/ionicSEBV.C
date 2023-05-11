
#include "ionicSEBV.h"

registerADMooseObject("ecmElectrochemApp", ionicSEBV);

InputParameters
ionicSEBV::validParams()
{

  InputParameters params = ADIntegratedBC::validParams();
  params.addClassDescription("Compute the outflow boundary condition.");

// Add a parameter with a default value; this value can be overridden in the input file.
  params.addParam<Real>("F_RT", 0.03868, "The constant of F/RT  when T = 300K, in unit 1/mV");
  params.addParam<Real>("reaction_rate", 0.5, "The reaction_rate of Li+ to Li, unitless.");
  params.addParam<Real>("LiPotElectrode", 0.0, "The Li potental in anode or cathode, in V.");
  params.addRequiredParam<Real>("ex_current", "The exchange current density at the interface, in unit mA/cm^2.");
  return params;
}

ionicSEBV::ionicSEBV(const InputParameters & parameters)
  : ADIntegratedBC(parameters),

  // Get the parameters from the input file
   _F_RT(getParam<Real>("F_RT")),
   _reaction_rate(getParam<Real>("reaction_rate")),
   _LiPotEle(getParam<Real>("LiPotElectrode")),
   _exchange_current(getParam<Real>("ex_current"))

   // Get the parameters from the input file
//    _ex_current(parameters.get<MaterialPropertyName>("ex_current")),
//    _exchange_current(getADMaterialProperty<Real>(_ex_current))
{
}

ADReal
ionicSEBV::computeQpResidual()
{
  ADReal k1 = std::exp(_reaction_rate * _F_RT * (_u[_qp] - _LiPotEle));
  ADReal k2 = std::exp(- (1 - _reaction_rate) * _F_RT * (_u[_qp] - _LiPotEle));
  return _test[_i][_qp] * _exchange_current * (k1 - k2);
}
