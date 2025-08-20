

#include "ADParticleSize.h"

registerMooseObject("ecmApp", ADParticleSize);
InputParameters
ADParticleSize::validParams()
{
    InputParameters params = ADMaterial::validParams();
    params.addClassDescription("Base class to define particle size or distribution");
    params.addParam<std::string>("base_name", "Base name for material class");
    params.addParam<MaterialPropertyName>("particle_size_name", "particle_size",
            "The name for the particle size material property");
    params.addParam<MaterialPropertyName>("surface_to_volume_name", "surface_to_volume",
            "The name for the surface_to_volume material property");
    params.addRequiredParam<Real>("particle_size", "Constant Particle size");
    params.addParam<FunctionName>("particle_size_function", "", "Function describing particle size");
    MooseEnum particleType("SPHERICAL CYLINDRICAL", "SPHERICAL");
    params.addParam<MooseEnum>("particleType", particleType, "Particle type for electrode");
    params.addParam<Real>("volume_fraction",1.0,"Volume fraction of particle in multi particle situations");
    return params;
}

ADParticleSize::ADParticleSize(const InputParameters & parameters)
        : ADMaterial(parameters),
        _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
        _particle_size(declareADPropertyByName<Real>(_base_name + getParam<MaterialPropertyName>("particle_size_name"))),
        _surface_to_volume(declareADPropertyByName<Real>(_base_name + getParam<MaterialPropertyName>("surface_to_volume_name"))),
        _size(getParam<Real>("particle_size")),
        _particle_size_function(getParam<FunctionName>("particle_size_function")!= "" ?
            &getFunction("particle_size_function") : NULL),
        _particle_type(getParam<MooseEnum>("particleType").getEnum<ParticleType>()),
        _volume_fraction(getParam<Real>("volume_fraction")),
        _vol_fraction(declareADPropertyByName<Real>(_base_name + "volume_fraction"))

{
    if (isParamValid("particle_size") && _particle_size_function)
    {
        mooseError("Cannot define both reaction rate and reaction rate function");
    }
    if (isParamValid("particle_size"))
    {
        if (_size <= 0)
            mooseError("particle size cannot be less than 0");
    }
}

void
ADParticleSize::computeQpProperties()
{
    if (_particle_size_function)
    {
        auto p = _q_point[_qp];
        _particle_size[_qp] = _particle_size_function->value(_t, p);
    } else
        _particle_size[_qp] = _size;
    _vol_fraction[_qp] = _volume_fraction;
    switch(_particle_type)
    {
        case ParticleType::Spherical:
            _surface_to_volume[_qp] = (3.0/_particle_size[_qp]);
            break;
        case ParticleType::Cylindrical:
            _surface_to_volume[_qp] = (2.0/_particle_size[_qp]);
            break;
    }
}

//Real
//ADParticleSize::surfacetoVolume()
//{
//    return MetaPhysicL::raw_value(3.0/_particle_size[_qp]);
//}
