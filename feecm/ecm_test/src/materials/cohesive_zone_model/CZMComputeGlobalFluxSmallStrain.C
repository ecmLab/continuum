

#include "CZMComputeGlobalFluxSmallStrain.h"


registerMooseObject("ecmApp", CZMComputeGlobalFluxSmallStrain);

InputParameters
CZMComputeGlobalFluxSmallStrain::validParams()
{
  InputParameters params = CZMcomputeGlobalFluxBase::validParams();

  params.addClassDescription(
      "Computes the czm flux in global coordinates for a small strain kinematic formulation");
  return params;
}

CZMComputeGlobalFluxSmallStrain::CZMComputeGlobalFluxSmallStrain(
    const InputParameters & parameters)
  : CZMcomputeGlobalFluxBase(parameters)
{
}

void
CZMComputeGlobalFluxSmallStrain::computeEquilibriumFluxAndDerivatives()
{
  _flux_global[_qp] = _interface_flux[_qp];
  _dflux_dvariablejump_global[_qp] = _dinterface_flux_dvariablejump[_qp];
  _dflux_djump_global[_qp] = _total_rotation[_qp] * _dinterface_flux_djump[_qp];

}
