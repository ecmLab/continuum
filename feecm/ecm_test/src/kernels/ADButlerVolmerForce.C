#include "ADButlerVolmerForce.h"

registerADMooseObject("ecmApp", ADButlerVolmerForce);

InputParameters
ADButlerVolmerForce::validParams()
{
    InputParameters params = ADKernel::validParams();
    params.addClassDescription("Implements a body force that is dependent on a material property");
    params.addParam<std::string>("base_name", "Base name for material class");
    params.addParam<Real>("scale", "Scaling coefficient");
    params.addParam<MaterialPropertyName>("mat_prop_name", "Name of material property");
    return params;
}

ADButlerVolmerForce::ADButlerVolmerForce(const InputParameters & parameters)
        : ADKernel(parameters),
        _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
        _scale(getParam<Real>("scale")),
        _mat_property(getADMaterialPropertyByName<Real>(_base_name + getParam<MaterialPropertyName>("mat_prop_name")))
{
}

ADReal
ADButlerVolmerForce::computeQpResidual()
{
    return _scale * _mat_property[_qp] * _test[_i][_qp];
}
