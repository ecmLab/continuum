

#include "IsotropicDiffusionMaterial.h"

registerADMooseObject("ecmApp", ADIsotropicDiffusionMaterial);

InputParameters ADIsotropicDiffusionMaterial::validParams()
{
    InputParameters params = Material::validParams();
    params.addClassDescription("Provides a diffusivity tensor based on 3 principal directions");
    params.addParam<std::string>("base_name", "Base name for material class");
    params.addRequiredParam<Real>("diffusion_coef", "Diffusivity");
    params.addParam<MaterialPropertyName>("diffusivity_name", "diffusivity",
            "Name of the diffusivity");
    return params;
}

ADIsotropicDiffusionMaterial::ADIsotropicDiffusionMaterial
    (const InputParameters & parameters) : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _diffusion_coef(getParam<Real>("diffusion_coef")),
    _diffusivity(declareADProperty<Real>(_base_name + getParam<MaterialPropertyName>("diffusivity_name")))
{

}

void ADIsotropicDiffusionMaterial::computeQpProperties()
{

    _diffusivity[_qp] = computeQpCorrection() * _diffusion_coef;
}
