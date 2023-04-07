
#include "SingleSEElectrodeBV.h"

registerADMooseObject("ecmElectrochemApp", SingleSEElectrodeBV);

InputParameters
SingleSEElectrodeBV::validParams()
{

  InputParameters params = ADIntegratedBC::validParams();
  params.addClassDescription("Compute the outflow boundary condition.");

// Add a parameter with a default value; this value can be overridden in the input file.
  params.addParam<Real>("F_RT", 38.68, "The constant of F/RT,in unit V, when T = 300K.");
  params.addParam<Real>("reaction_rate", 0.5, "The reaction_rate of Li+ to Li, unitless.");
  params.addParam<Real>("LiPotElectrode", 0.0, "The Li potental in anode or cathode, in V.");
  params.addRequiredParam<Real>("exchange_current", "The exchange current density at the interface ($\\mathrm{A/cm^2}$).");
  return params;
}

SingleSEElectrodeBV::SingleSEElectrodeBV(const InputParameters & parameters)
  : ADIntegratedBC(parameters),

  // Get the parameters from the input file
   _F_RT(getParam<Real>("F_RT")),
   _reaction_rate(getParam<Real>("reaction_rate")),
   _LiPotEle(getParam<Real>("LiPotElectrode")),
   _exchange_current(getParam<Real>("exchange_current"))
{
}

ADReal
SingleSEElectrodeBV::computeQpResidual()
{
  ADReal k1 = std::exp(_reaction_rate * _F_RT * (_u[_qp] - _LiPotEle));
  ADReal k2 = std::exp(- (1 - _reaction_rate) * _F_RT * (_u[_qp] - _LiPotEle));
  return _test[_i][_qp] * _exchange_current * (k1 - k2);
}
