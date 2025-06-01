

#include "ADCoupledButlerVolmerForce.h"

registerMooseObject("ecmApp",ADCoupledButlerVolmerForce);

InputParameters
ADCoupledButlerVolmerForce::validParams()
{
    InputParameters params = ADKernel::validParams();
    params.addClassDescription("Implements a butler volmer type source term");
    params.addRequiredCoupledVar("electrolyte_potential", "Variable for "
                    "potenital of the electrolyte");
    params.addRequiredCoupledVar("electrode_potential", "Variable for p"
                    "potenital of the electrode");
    params.addParam<Real>("faraday",96.4853329, "Faraday's constant");
    params.addParam<Real>("temperature", 298, "Temperature");
    params.addParam<Real>("gas_constant", 8.314462681,"Universal Gas Constant");
    params.addParam<Real>("exchange_current_density",
            "Constant exchange current density");
    return params;
}

ADCoupledButlerVolmerForce::ADCoupledButlerVolmerForce(const InputParameters & parameters)
       : ADKernel(parameters),
        _i0(getParam<Real>("exchange_current_density")),
        _faraday(getParam<Real>("faraday")),
        _gas_constant(getParam<Real>("gas_constant")),
        _temperature(getParam<Real>("temperature")),
        _num_electrolyte_potential_var(coupled("electrolyte_potential_var")),
        _electrolyte_potential(adCoupledValue("electrolyte_potential_var")),
        _num_electrode_potential_var(coupled("electrode_potential_var")),
        _electrode_potential(adCoupledValue("electrode_potential_var")),
        _equilibrium_potential(getADMaterialProperty<Real>("equlibrium_potential"))
{
}

ADReal
ADCoupledButlerVolmerForce::computeQpResidual()
{
    auto RTF = _gas_constant * _temperature / _faraday;
    auto eta = _electrode_potential[_qp] - _electrolyte_potential[_qp] -
                _equilibrium_potential[_qp];
    return 2.0 * _i0 * std::sinh(0.5*eta/RTF);
}
