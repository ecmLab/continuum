
#include "MixSEIntBV.h"

registerADMooseObject("ecmApp", MixSEIntBV);

InputParameters
MixSEIntBV::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addClassDescription("Compute the outflow boundary condition.");
  params.addRequiredParam<MaterialPropertyName>("LiPotElectrode", "The Li potental in anode or cathode, in V.");

// Add a coupled parameter: potEn
  params.addRequiredCoupledVar("potEn", "The potential of Electron");

// Add a parameter with a default value; this value can be overridden in the input file.
    params.addParam<Real>(
        "F_RT",
        38.68,
        "The constant of F/RT,in unit V, when T = 300K.");

  return params;
}

MixSEIntBV::MixSEIntBV(const InputParameters & parameters)
  : ADIntegratedBC(parameters),
  // Couple to the potential of electron
   _potEn(adCoupledValue("potEn")),

  // Get the parameters from the Material object
   _exchange_current(getADMaterialProperty<Real>("exchange_current")),
   _reaction_rate(getADMaterialProperty<Real>("reaction_rate")),

  // Get the parameters from the input file
   _LiPotElectrode(parameters.get<MaterialPropertyName>("LiPotElectrode")),
   _LiPotEle(getADMaterialProperty<Real>(_LiPotElectrode)),
   _F_RT(getParam<Real>("F_RT"))
{
}

ADReal
MixSEIntBV::computeQpResidual()
{
  ADReal k1 = std::exp(_reaction_rate[_qp] * _F_RT * (_u[_qp] + _potEn[_qp] - _LiPotEle[_qp]));
  ADReal k2 = std::exp(- (1 - _reaction_rate[_qp]) * _F_RT * (_u[_qp] + _potEn[_qp] - _LiPotEle[_qp]));
  return _test[_i][_qp] * _exchange_current[_qp] * (k1 - k2);
}
