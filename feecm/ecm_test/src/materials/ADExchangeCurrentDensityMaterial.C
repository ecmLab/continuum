

#include "ADExchangeCurrentDensityMaterial.h"

registerADMooseObject("ecmApp", ADExchangeCurrentDensityMaterial);

InputParameters
ADExchangeCurrentDensityMaterial::validParams()
{
    InputParameters params = ADMaterial::validParams();
    params.addClassDescription("Material that computes the exchange current density");
    params.addParam<std::string>("base_name", "Base name for material class");
    MooseEnum currentType("CONSTANT LI_INSERTION", "CONSTANT");
    params.addParam<MooseEnum>("exchange_current_density_type",currentType, "Exchange current density type" );
    params.addRequiredCoupledVar("concentration", "Concentration variable");
    params.addParam<MaterialPropertyName>("exchange_current_density_name",
            "exchange_current_density", "Name of the exchange current density");
    params.addParam<Real>("reference_current_density", "Constant Exchange current density");
    params.addParam<Real>("cmax", "Maximum concentration");
    return params;
}

ADExchangeCurrentDensityMaterial::ADExchangeCurrentDensityMaterial(const InputParameters & parameters)
        : ADMaterial(parameters),
        _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
        _i0(declareADProperty<Real>(_base_name + getParam<MaterialPropertyName>("exchange_current_density_name"))),
        _exchange_type(getParam<MooseEnum>("exchange_current_density_type").getEnum<ExchangeCurrentDensityType>()),
        _concentration(adCoupledValue("concentration")),
        _i_ref(getParam<Real>("reference_current_density")),
        _cmax(getParam<Real>("cmax"))
{
    if (_exchange_type == ExchangeCurrentDensityType::Lithium_Insertion && _cmax <= 0.0)
        mooseError("Need to have a positive concentration value when using Li insertion kinetics");
}

void
ADExchangeCurrentDensityMaterial::computeQpProperties()
{
    switch(_exchange_type)
    {
        case ExchangeCurrentDensityType::Constant:
            _i0[_qp] = _i_ref;
            break;
        case ExchangeCurrentDensityType::Lithium_Insertion:
            auto cbar = _concentration[_qp] /_cmax;
            if (_concentration[_qp] < 0.001 * _cmax)
                cbar = 0.001;
            if (_concentration[_qp] > 0.999 * _cmax)
                cbar = 0.999;
            _i0[_qp] = 2.0 * _i_ref * std::sqrt(cbar * (1.0 - cbar));
            break;
    }
}
