
#include "MixSEMetalBV.h"

registerADMooseObject("ecmApp", MixSEMetalBV);

InputParameters
MixSEMetalBV::validParams()
{

   InputParameters params = ADIntegratedBC::validParams();
   params.addClassDescription("Compute the deposition boundary condition.");

// Add a coupled parameter: potLi
   params.addRequiredCoupledVar("potRef", "The potential of the other species");

// Add a parameter with a default value; this value can be overridden in the input file.
   params.addParam<Real>(
        "F_RT",
        38.68,
        "The constant of F/RT, when T = 300K.");

   return params;
}

MixSEMetalBV::MixSEMetalBV(const InputParameters & parameters)
  : ADIntegratedBC(parameters),
  // Couple to the other potential
   _potRef(adCoupledValue("potRef")),

  // Get the parameters from the input file
   _exchange_current(getADMaterialProperty<Real>("exchange_current")),
   _reaction_rate(getADMaterialProperty<Real>("reaction_rate")),
   _electron_concentration(getADMaterialProperty<Real>("electron_concentration")),
   _F_RT(getParam<Real>("F_RT"))
{
}

ADReal
MixSEMetalBV::computeQpResidual()
{
  ADReal k1 = std::exp(_reaction_rate[_qp] * _F_RT * (_potRef[_qp] + _u[_qp] ));
  ADReal k2 = std::exp(- (1 - _reaction_rate[_qp]) * _F_RT * (_potRef[_qp] + _u[_qp] ));
//  ADReal k1 = _electron_concentration[_qp] * _exchange_current[_qp] * _F_RT * (_potRef[_qp] + _u[_qp]);

  return _test[_i][_qp] * _electron_concentration[_qp] * _exchange_current[_qp] * (k1 - k2);
//  return _test[_i][_qp] * k1;
}
