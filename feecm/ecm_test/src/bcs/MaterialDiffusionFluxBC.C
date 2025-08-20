
#include "MaterialDiffusionFluxBC.h"


registerMooseObject("ecmApp", MaterialDiffusionFluxBC);

InputParameters
MaterialDiffusionFluxBC::validParams()
{
    InputParameters params = ADIntegratedBC::validParams();
    params.addParam<MaterialPropertyName>("diffusivity", "diffusivity", "Name of diffusivity");
    return params;
}


MaterialDiffusionFluxBC::MaterialDiffusionFluxBC(const InputParameters& parameters)
        : ADIntegratedBC(parameters),
        _diffusivity(&getADMaterialProperty<Real>("diffusivity"))
{

}


ADReal
MaterialDiffusionFluxBC::computeQpResidual()
{
    return -(*_diffusivity)[_qp] * _grad_u[_qp] * _normals[_qp] * _test[_i][_qp];
}
