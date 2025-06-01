

#include "ComputeSurfaceConcentrationAux.h"

registerMooseObject("ecmApp", ComputeSurfaceConcentrationAux);

InputParameters
ComputeSurfaceConcentrationAux::validParams()
{
    InputParameters params = AuxKernel::validParams();
    params.addClassDescription("Compute surface concentration from average concentration"
            "based on 2 term approximation");
    params.addRequiredCoupledVar("concentration", "Variable representing avg concentration");
    params.addRequiredParam<MaterialPropertyName>("flux", "Name of the flux variable");
    params.addRequiredParam<Real>("particle_size","Size of particle");
    params.addRequiredParam<MaterialPropertyName>("diffusion_coefficient", "diffusion coefficient");
    params.addRequiredParam<Real>("scale", "scale");
    return params;
}

ComputeSurfaceConcentrationAux::ComputeSurfaceConcentrationAux(const InputParameters & parameters)
    : AuxKernel(parameters),
        _concentration(adCoupledValue("concentration")),
        _flux(getADMaterialProperty<Real>("flux")),
        _particle_size(getParam<Real>("particle_size")),
        _diffusion_coefficient(getADMaterialProperty<Real>("diffusion_coefficient")),
        _scale(getParam<Real>("scale"))
{

}

Real
ComputeSurfaceConcentrationAux::computeValue()
{
    return MetaPhysicL::raw_value(_concentration[_qp] + _particle_size * _flux[_qp]/_scale/_diffusion_coefficient[_qp]);
//    return MetaPhysicL::raw_value(_concentration[_qp]);
}
