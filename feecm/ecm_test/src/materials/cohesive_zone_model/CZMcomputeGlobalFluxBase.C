
#include "CZMcomputeGlobalFluxBase.h"


InputParameters
CZMcomputeGlobalFluxBase::validParams()
{
  InputParameters params = InterfaceMaterial::validParams();

  params.addClassDescription(
      "Base class for computing the equilibrium flux and its derivatives.");
  params.suppressParameter<bool>("use_displaced_mesh");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addParam<bool>("include_gap", true, "include gap displ");
  return params;
}

CZMcomputeGlobalFluxBase::CZMcomputeGlobalFluxBase(const InputParameters & parameters)
  : InterfaceMaterial(parameters),
    _base_name(isParamValid("base_name") && !getParam<std::string>("base_name").empty()
                   ? getParam<std::string>("base_name") + "_"
                   : ""),
    _include_gap(getParam<bool>("include_gap")),
    _flux_global(declarePropertyByName<Real>(_base_name + "flux_global")),
    _interface_flux(
        getMaterialPropertyByName<Real>(_base_name + "interface_flux")),
    _dflux_dvariablejump_global(
        declarePropertyByName<Real>(_base_name + "dflux_dvariablejump_global")),
    _dflux_djump_global(
        declarePropertyByName<RealVectorValue>(_base_name + "dflux_djump_global")),
    _dinterface_flux_dvariablejump(
        getMaterialPropertyByName<Real>(_base_name + "dinterface_flux_dvariablejump")),
    _dinterface_flux_djump(
        getMaterialPropertyByName<RealVectorValue>(_base_name + "dinterface_flux_djump")),
    _total_rotation(getMaterialPropertyByName<RankTwoTensor>(_base_name + "total_rotation"))

{
}

void
CZMcomputeGlobalFluxBase::computeQpProperties()
{
  // rotate local flux and derivatives to the global coordinate system
  computeEquilibriumFluxAndDerivatives();
}
