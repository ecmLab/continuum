

#include "ADComputeEquilibriumPotential.h"
#include "Function.h"

registerADMooseObject("ecmApp", ADComputeEquilibriumPotential);

InputParameters
ADComputeEquilibriumPotential::validParams()
{
    InputParameters params = ADMaterial::validParams();
    params.addClassDescription(" Computes the equilibrium potential and chemical potential"
                    "from the swelling_chemical_potential,"
                    "elastic_chemical_potential, concentration and reaction rate. ");
    params.addParam<std::string>("base_name", "Base name for material class");
    params.addParam<Real>("faraday",96.4853329, "Faraday's constant");
    params.addParam<Real>("temperature", 298, "Temperature");
    params.addParam<Real>("R", 8.314462681,"Universal Gas Constant");
    params.addParam<bool>("include_reaction_rate", true, "Include reaction rate"
                           "Nearnst potential to equilibrium potential");
    params.addParam<bool>("include_conc", false, "Include concentration term to Nearnst potential");
    params.addParam<Real>("reaction_rate", "Reaction rate ");
    params.addParam<FunctionName>("reaction_rate_function", "", "Function describing reaction rate");
    params.addParam<bool>("is_ocv",true, "Is the function actually an OCV curve");
    params.addParam<bool>("include_mechanical_effects", false, "Include mechanics contributions to OCV");
    params.addParam<bool>("exclude_elastic_contribution", true, "Exclude elastic energy contributions to OCV");
    params.addCoupledVar("concentration", "Coupled concentration");
    params.addParam<Real>("cref", "Reference Concentration");
    MooseEnum potentialType("DEHOFF BOWER", "BOWER");
    params.addParam<MooseEnum>("potential", potentialType, "Type of diffusion potential");
    params.addParam<Real>("nucleation_overpotential", 0.0, "Nucleation overpotential for this material");
    return params;
}

ADComputeEquilibriumPotential::ADComputeEquilibriumPotential(const InputParameters & parameters)
    : ADMaterial(parameters),
      _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
      _faraday(getParam<Real>("faraday")),
      _temp(getParam<Real>("temperature")),
      _gas_constant(getParam<Real>("R")),
      _include_reaction_rate(getParam<bool>("include_reaction_rate")),
      _has_conc(isCoupled("concentration")),
      _include_conc(getParam<bool>("include_conc")),
      _concentration(isParamValid("concentration") ?  &adCoupledValue("concentration") : nullptr),
      _c0(getParam<Real>("cref")),
      _include_mechanical_effects(getParam<bool>("include_mechanical_effects")),
      _exclude_elastic_energy(getParam<bool>("exclude_elastic_contribution")),
      _reaction_rate(isParamValid("reaction_rate") ? getParam<Real>("reaction_rate") : 0),
      _reaction_rate_function(getParam<FunctionName>("reaction_rate_function")!= ""
            ? &getFunction("reaction_rate_function") : NULL),
      _is_ocv(getParam<bool>("is_ocv")),
      _swelling_chemical_potential(_include_mechanical_effects ?
          &getADMaterialProperty<Real>(_base_name + "swelling_chemical_potential")
           : nullptr),
      _elastic_chemical_potential(_include_mechanical_effects ?
          &getADMaterialProperty<Real>(_base_name + "elastic_chemical_potential")
        : nullptr),
      _equilibrium_potential(declareADPropertyByName<Real>(_base_name + "equilibrium_potential")),
      _chemical_potential(declareADProperty<Real>(_base_name + "chemical_potential")),
      _state_of_charge(declareADProperty<Real>(_base_name + "state_of_charge")),
      _diff_potential(getParam<MooseEnum>("potential").getEnum<DiffPotential>())
{

    if (isParamValid("reaction_rate") && _reaction_rate_function)
    {
        mooseError("Cannot define both reaction rate and reaction rate function");
    }
}

void
ADComputeEquilibriumPotential::computeQpProperties()
{
    auto prefac = _gas_constant * _temp;
    auto prefac2 = 1.0/_faraday;
    _equilibrium_potential[_qp] = 0.0;
    _chemical_potential[_qp] = 0.0;
    ADReal cbar = 0.0;
    ADReal c = 0.0;
    if (_has_conc) {
        c = (*_concentration)[_qp] > 0 ? (*_concentration)[_qp]: ADReal(0.0);
        cbar = c/_c0;
        if ((*_concentration)[_qp] < 0.001 * _c0)
                cbar = 0.001;
        if ((*_concentration)[_qp] > 0.999 * _c0)
            cbar = 0.999;
    }
    _state_of_charge[_qp]  = cbar;


    if (_include_reaction_rate)
    {
        if (_reaction_rate_function)
        {
//            auto p = _q_point[_qp];
            const Point p;
            auto reaction_rate = _reaction_rate_function->value(MetaPhysicL::raw_value(cbar), p);
            _equilibrium_potential[_qp] += (_is_ocv ? reaction_rate : reaction_rate*prefac*prefac2);
        } else
        {
            _equilibrium_potential[_qp] += (_is_ocv ? _reaction_rate : prefac*prefac2*_reaction_rate);
        }
    }

    if (_has_conc && _include_conc)
    {
        ADReal chem = 0.0;
        if ((*_concentration)[_qp] > 0)
        {
            if (_diff_potential == DiffPotential::DeHoff)
            {
//            {   if ((*_concentration)[_qp] > 0 && (_c0 - (*_concentration)[_qp] > 0))
                if (c > 0 && (_c0 - c) > 0)
                    chem = std::log(_c0/(_c0 - c));
            }
            else if (_diff_potential == DiffPotential::Bower)
                chem = std::log(cbar);

//            _equilibrium_potential[_qp] -= prefac * chem * prefac2;
            _chemical_potential[_qp] += prefac * chem;
        }
    }
    if (_include_mechanical_effects)
    {
        auto mech = (*_swelling_chemical_potential)[_qp];
        if (!_exclude_elastic_energy)
            mech += (*_elastic_chemical_potential)[_qp];
//        _equilibrium_potential[_qp] -=  mech * prefac2;
        _chemical_potential[_qp] += mech;
    }
    _equilibrium_potential[_qp] -= prefac2 * _chemical_potential[_qp];

}
