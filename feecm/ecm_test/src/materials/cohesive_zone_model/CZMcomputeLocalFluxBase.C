
#include "Assembly.h"
#include "CZMcomputeLocalFluxBase.h"

InputParameters
CZMcomputeLocalFluxBase::validParams()
{
  InputParameters params = InterfaceMaterial::validParams();

  params.addClassDescription("Base class for implementing cohesive zone constitutive material "
                             "models that can be formulated using the total displacement jump");
  params.addRequiredCoupledVar("variable",
                               "The string of variable suitable for the problem statement");
  params.suppressParameter<bool>("use_displaced_mesh");
  params.addParam<std::string>("base_name", "Material property base name");
      params.addParam<bool>("include_gap", true, "include displacement"
        "contribution to the off diagonal jacobian, setting to false will only "
        "include diagonal jacobian components");
  return params;
}

CZMcomputeLocalFluxBase::CZMcomputeLocalFluxBase(const InputParameters & parameters)
  : InterfaceMaterial(parameters),
    _base_name(isParamValid("base_name") && !getParam<std::string>("base_name").empty()
                   ? getParam<std::string>("base_name") + "_"
                   : ""),
    _include_gap(getParam<bool>("include_gap")),
    _interface_flux(declarePropertyByName<Real>(_base_name + "interface_flux")),
     _dinterface_flux_dvariablejump(
        declarePropertyByName<Real>(_base_name + "dinterface_flux_dvariablejump")),
    _dinterface_flux_djump(
        declarePropertyByName<RealVectorValue>(_base_name + "dinterface_flux_djump")),
    _interface_variable_jump(
        getMaterialPropertyByName<Real>(_base_name + "interface_variable_jump")),
    _interface_displacement_jump(_include_gap ?
        &getMaterialPropertyByName<RealVectorValue>(_base_name + "interface_displacement_jump") : nullptr)
{
}

void
CZMcomputeLocalFluxBase::initQpStatefulProperties()
{
  _interface_flux[_qp] = 0;
}

void
CZMcomputeLocalFluxBase::computeQpProperties()
{
  computeInterfaceFluxAndDerivatives();
}
