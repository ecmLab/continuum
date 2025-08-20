

#include "CZMButlerVolmer.h"
#include "Function.h"

registerMooseObject("ecmApp", CZMButlerVolmer);

InputParameters
CZMButlerVolmer::validParams()
{
    InputParameters params = CZMcomputeLocalFluxBase::validParams();
    params.addClassDescription("Material class that computes a flux based on "
        "jump of diffusive type variable. The class also computes the residual "
        "and jacobian contribution based on simple linear model of flux with "
        "with interface separation flux = (flux * (d -interface_jump(0)");
    params.addRequiredParam<FunctionName>("k_function", "Function of the relevant gap problem");
    params.addRequiredParam<Real>("max_separation","separation at failure for 0 flux");
    MooseEnum conductanceType("CONDUCTANCE RCT EXCHANGE_CURRENT_DENSITY", "EXCHANGE_CURRENT_DENSITY");
    params.addParam<MooseEnum>("conductanceType", conductanceType, "What type of conductance");
    params.addCoupledVar("eq_potential_coupled_var", "The name of coupled variable "
        "controlling the equilibrium potential");
    params.addParam<bool>("include_equil", false, "Include contribution of "
    "the equilibrium potential");
    MooseEnum surfaceType("PRIMARY SECONDARY", "SECONDARY");
    params.addParam<MooseEnum>("surfaceType", surfaceType, "Electrode surface");
    params.addParam<Real>("cref",1.0,"Reference Concentration for equilibrium potential");
    params.addParam<Real>("faraday",96.4853329, "Faraday's constant");
    params.addParam<Real>("temperature", 298, "Temperature");
    params.addParam<Real>("R", 8.314462681,"Universal Gas Constant");
    params.addParam<FunctionName>("reaction_rate_function", "", "Function describing reaction rate");

    MooseEnum computeType("LINEAR_BUTLER_VOLMER BUTLER_VOLMER LITHIUM_INSERTION", "LINEAR_BUTLER_VOLMER");
    params.addParam<MooseEnum>("computeType", computeType, "What do we want this constraint to do");

    return params;
}

CZMButlerVolmer::CZMButlerVolmer(const InputParameters& parameters)
        : CZMcomputeLocalFluxBase(parameters),
        _d(getParam<Real>("max_separation")),
        _k_function(&this->getFunction("k_function")),
        _conductance(getParam<MooseEnum>("conductanceType").getEnum<Conductance>()),
        _include_equil(getParam<bool>("include_equil")),
        _surface_type(getParam<MooseEnum>("surfaceType").getEnum<SurfaceType>()),
        _compute_type(getParam<MooseEnum>("computeType").getEnum<Compute>()),
        _cref(getParam<Real>("cref")),
        _reaction_rate_function(getParam<FunctionName>("reaction_rate_function")!= ""
            ? &getFunction("reaction_rate_function") : NULL),
        _faraday(getParam<Real>("faraday")),
        _temperature(getParam<Real>("temperature")),
        _gas_constant(getParam<Real>("R")),
        _equilibrium_potential(declarePropertyByName<Real>("equilibrium_potential")),
        _RTF(_gas_constant * _temperature/_faraday),
        _flux_sign(_surface_type == SurfaceType::Secondary ? 1.0 : -1.0)

{

    if (_include_equil)
    {
        if (_surface_type == SurfaceType::Secondary)
        {
            _coupled_value = &coupledNeighborValue("eq_potential_coupled_var",0);
        }
        else if (_surface_type == SurfaceType::Primary) {
            _coupled_value = &coupledValue("eq_potential_coupled_var",0);
        } else
        {
            _coupled_value = nullptr;
        }
        if (!(_coupled_value))
            mooseError("Coupled variable must be given to compute equilibrium potential");
    } else {
        _coupled_value = nullptr;
        if (_coupled_value)
            mooseWarning("coupled var is provided it is not used in anything. Perhaps provide reaction rate function");
    }

}

void
CZMButlerVolmer::computeInterfaceFluxAndDerivatives()
{
    // interface jump is always computed as secondary - primary
    // In the butler volmer problem the overpotential is always defined as Electrode - Electrolyte
    // So in the material class we can adjust which surface is the electrode, using surface_type parameter
    /// The default for this parameter
    Real cbar = 0.0;

    _equilibrium_potential[_qp] = 0.0;
    if (_coupled_value){
        cbar = (*_coupled_value)[_qp]/_cref;
        const Point p;
        _equilibrium_potential[_qp] = _reaction_rate_function->value(cbar, p);
    }
    // Compute value of gap conductance from the function
    // TODO can make this a function of space
    const Point p;
    auto k = _k_function->value(_interface_variable_jump[_qp], p);
    auto i0 = k;
    switch(_conductance)
    {
        /// Convert all quantities to exchange current density for uniform formulation
        case Conductance::Rct:
            if (k > 0.0)
                i0 = _RTF /k;
            else i0 = 1e20;
            break;
        case Conductance::Conductance:
            i0 *= _RTF;
        case Conductance::Exchange_Current_density:
            i0 = k;

    }

    auto flux = i0;
    auto flux_deriv = i0;
    if (_compute_type == Compute::Linear_Butler_Volmer)
    {
        flux *= (_flux_sign * _interface_variable_jump[_qp]
                       - _equilibrium_potential[_qp])/_RTF;
        flux_deriv = i0 * _flux_sign / _RTF;
    }
    else if (_compute_type == Compute::Butler_Volmer)
    {
        flux = i0 * 2.0 * std::sinh(( _flux_sign * _interface_variable_jump[_qp]
                                        - _equilibrium_potential[_qp])/_RTF/2.0);
        flux_deriv = i0 * 2.0 * _flux_sign * std::cosh((_flux_sign * _interface_variable_jump[_qp]
                                        - _equilibrium_potential[_qp])/_RTF/2.0)
                                        /_RTF/2.0;
    }
    else {
        flux = 0.0;
        flux_deriv = 0.0;
    }

    if (_include_gap)
    {
        // Include other models here
        if ((*_interface_displacement_jump)[_qp](0) < _d){
           _interface_flux[_qp] = flux * (-(*_interface_displacement_jump)[_qp](0) + _d );
            _dinterface_flux_dvariablejump[_qp] = flux_deriv * (-(*_interface_displacement_jump)[_qp](0) + _d );

            _dinterface_flux_djump[_qp](0) = -flux;
            _dinterface_flux_djump[_qp](1) = 0;
            _dinterface_flux_djump[_qp](2) = 0;
        } else
        {
            _interface_flux[_qp] = 0.0;
            _dinterface_flux_dvariablejump[_qp] = 0.0;
            _dinterface_flux_djump[_qp](0) = 0.0;
            _dinterface_flux_djump[_qp](1) = 0;
            _dinterface_flux_djump[_qp](2) = 0;

        }
    } else
    {
        _interface_flux[_qp] = flux;
        _dinterface_flux_dvariablejump[_qp] = flux_deriv;
        _dinterface_flux_djump[_qp](0) = 0.0;
        _dinterface_flux_djump[_qp](1) = 0;
        _dinterface_flux_djump[_qp](2) = 0;
    }

}
