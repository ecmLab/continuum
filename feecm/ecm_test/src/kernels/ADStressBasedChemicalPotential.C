

#include "ADStressBasedChemicalPotential.h"

registerMooseObject("ecmApp", ADStressBasedChemicalPotential);

InputParameters
ADStressBasedChemicalPotential::validParams()
{
    InputParameters params = ADKernel::validParams();
    params.addClassDescription("Stress based chemical potential, calculates the "
            "the contribution of elastic energy and swelling based energy to the "
            "chemical potential. These quantities are all calculated in the reference"
            "space");
    params.addParam<std::string>("base_name", "Material property base name");
    params.set<bool>("use_displaced_mesh", false);
    return params;
}

ADStressBasedChemicalPotential::ADStressBasedChemicalPotential(const InputParameters & parameters)
    : ADKernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _swelling_chemical_potential(getADMaterialProperty<Real>(_base_name + "swelling_chemical_potential")),
    _elastic_chemical_potential(getADMaterialProperty<Real>(_base_name + "elastic_chemical_potential"))
{

}

ADReal
ADStressBasedChemicalPotential::computeQpResidual()
{
    return (_u[_qp] - _swelling_chemical_potential[_qp] -
            _elastic_chemical_potential[_qp]) * _test[_i][_qp];
}
