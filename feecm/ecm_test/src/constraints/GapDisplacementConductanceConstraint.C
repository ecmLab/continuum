

#include "GapDisplacementConductanceConstraint.h"
#include "Function.h"
#include "SubProblem.h"
#include "MaterialBase.h"

registerADMooseObject("ecmApp", GapDisplacementConductanceConstraint);

InputParameters
GapDisplacementConductanceConstraint::validParams()
{
	InputParameters params = ADMortarConstraint::validParams();
    params.addClassDescription(
        "Computes the residual and Jacobian contributions for the 'Lagrange Multiplier' "
        "implementation of the thermal contact problem. For more information, see the "
        "detailed description here: http://tinyurl.com/gmmhbe9");
    params.addParam<FunctionName>("k_function", "Function of the relevant gap problem");
    MooseEnum conductanceType("CONDUCTANCE RCT EXCHANGE_CURRENT_DENSITY", "CONDUCTANCE");
    params.addParam<MooseEnum>("conductanceType", conductanceType, "What type of conductance");
    params.addParam<Real>("k", "Gap conductance");
    params.addParam<Real>("faraday", 96485.3329, "Faraday's Constant");
    params.addParam<Real>("R", 8.3145, "Universal Gas Constant");
    params.addParam<Real>("temperature", 298, "Value of temperature to use");
    params.addCoupledVar("displacements", "Displacement variables");
    params.addParam<NonlinearVariableName>("concentration", "Concentration variable name");
    params.addParam<bool>("include_equilibrium_potential", false, "Equilibrium Potential of Electrode");
    params.addParam<bool>("include_gap", false, "Include gap value in calculation");
    params.addParam<bool>("include_concentration", true, "Include concentration value in calculation");
    MooseEnum surfaceType("PRIMARY SECONDARY", "SECONDARY");
    params.addParam<MooseEnum>("surfaceType", surfaceType, "Electrode surface");
    MooseEnum computeType("LINEAR_BUTLER_VOLMER BUTLER_VOLMER", "LINEAR_BUTLER_VOLMER");
    params.addParam<MooseEnum>("computeType", computeType, "What do we want this constraint to do");
    return params;
}


GapDisplacementConductanceConstraint::GapDisplacementConductanceConstraint(const InputParameters & parameters)
  : ADMortarConstraint(parameters),
    _k_function(isParamValid("k_function") ? & this->getFunction("k_function") : NULL),
    _conductance(getParam<MooseEnum>("conductanceType").getEnum<Conductance>()),
    _my_k(getParam<Real>("k")),
    _gas_constant(getParam<Real>("R")),
    _temp(getParam<Real>("temperature")),
    _faraday(getParam<Real>("faraday")),
    _disp_name(parameters.getVecMooseType("displacements")),
    _n_disp(_disp_name.size()),
    _disp_secondary(_n_disp),
    _disp_primary(_n_disp),
    _include_concentration(getParam<bool>("include_concentration")),
    _concentration_var((isParamValid("concentration") && _include_concentration)

                          ? &this->_subproblem.getStandardVariable(
                                _tid, parameters.getMooseType("concentration"))
                          : nullptr),
    _include_equilibrium_potential(getParam<bool>("include_equilibrium_potential")),
    _include_gap(getParam<bool>("include_gap")),
    _surface_type(getParam<MooseEnum>("surfaceType").getEnum<SurfaceType>()),
    _compute(getParam<MooseEnum>("computeType").getEnum<Compute>()),
    _RTF(_gas_constant * _temp/_faraday)

{
    if (_k_function == NULL && !parameters.isParamSetByUser("k"))
        mooseError("Either k or k_function must be given");
    if (_surface_type == SurfaceType::Primary)
    {
        _flux_sign = -1.0;
        if (_concentration_var)  _concentration = &_concentration_var->adSlnNeighbor();

    }
    else if (_surface_type == SurfaceType::Secondary)
    {
        _flux_sign = 1.0;
        if (_concentration_var)  _concentration = &_concentration_var->adSln();
    }
//    if (_concentration_var)
//    {
//        /// (TBD) Need to add code to make sure that primary and secondary surfaces can
//        //// be switched
//        if (_surface_type == SurfaceType::Primary)
//        {
//            _concentration = &_concentration_var->adSlnNeighbor();
//            _flux_sign = -1.0;
//        }
//        else if (_surface_type == SurfaceType::Secondary)
//        {
//            _concentration = &_concentration_var->adSln();
//            _flux_sign = 1.0;
//        }
//    }
    if (_include_gap)
    {
        for (unsigned int i = 0; i < _n_disp; ++i)
        {
            auto disp_var = &this->_subproblem.getStandardVariable(_tid, _disp_name[i]);
            _disp_secondary[i] = &disp_var->adSln();
            _disp_primary[i] = &disp_var->adSlnNeighbor();
        }
    }
    if (_include_equilibrium_potential)
    {
        if (_surface_type == SurfaceType::Primary)
        {
            _equilibrium_potential = (&getNeighborADMaterialProperty<Real>("equilibrium_potential") ?
                                        &getNeighborADMaterialProperty<Real>("equilibrium_potential") : nullptr);
        }
        else if (_surface_type == SurfaceType::Secondary) {
            _equilibrium_potential = (&getADMaterialProperty<Real>("equilibrium_potential") ?
                                        &getADMaterialProperty<Real>("equilibrium_potential") : nullptr);
        }
    }
}

ADReal
GapDisplacementConductanceConstraint::computeQpResidual(Moose::MortarType mortar_type)
{

  switch (mortar_type)
  {
    case Moose::MortarType::Primary:
      return _lambda[_qp] * _test_primary[_i][_qp];
    case Moose::MortarType::Secondary:
      return -_lambda[_qp] * _test_secondary[_i][_qp];
    case Moose::MortarType::Lower:
    {
        auto residual = _lambda[_qp] * _test[_i][_qp];
        _k = _my_k;
        auto i0 = _k;
//        if (_has_primary)
        {
            if (_include_gap)
            {
                // we are creating a dual version of phys points primary and secondary here...
                ADRealVectorValue dual_phys_points_primary;
                ADRealVectorValue dual_phys_points_secondary;
                ADRealVectorValue dual_normals;
                for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
                {
                  dual_phys_points_primary(i) = _phys_points_primary[_qp](i);
                  dual_phys_points_secondary(i) = _phys_points_secondary[_qp](i);
                  dual_normals(i) = _normals[_qp](i);
                }

                // ...which uses the derivative vector of the primary and secondary displacements as
                // an approximation of the true phys points derivatives
                for (unsigned int i = 0; i < _n_disp; ++i)
                {
                  dual_phys_points_primary(i).derivatives() = (*_disp_primary[i])[_qp].derivatives();
                  dual_phys_points_secondary(i).derivatives() = (*_disp_secondary[i])[_qp].derivatives();
                }

                auto gap = (dual_phys_points_primary - dual_phys_points_secondary) * dual_normals;
                if (_k_function)
                    computeConductance(gap);
            }

            switch(_conductance)
            {
                /// Convert all quantities to exchange current density for uniform formulation
                case Conductance::Rct:
                    if (_k > 0.0)
                        i0 = _RTF /_k;
                    else i0 = 1e20;
                    break;
                case Conductance::Conductance:
                    i0 *= _RTF;
            }
            bool check_conc = (_concentration && _include_concentration && (*_concentration)[_qp] < 0);
            auto eq_pot_qp = ((_include_equilibrium_potential && _equilibrium_potential)
                            ? (*_equilibrium_potential)[_qp] : ADReal(0.0));
            if (check_conc)
                return residual;
            else
            {
                /// (TBD) check to make sure sign is right when primary and secondary
                /// surfaces are switched ....
                if (_compute == Compute::Linear_Butler_Volmer)
                    residual = _test[_i][_qp] *
                        (_lambda[_qp] +  i0/_RTF * (_flux_sign * (_u_secondary[_qp] - _u_primary[_qp]) - eq_pot_qp));
                else if (_compute == Compute::Butler_Volmer)
                {
                    residual = _test[_i][_qp] *
                        (_lambda[_qp] + 2.0 *i0 *std::sinh(( _flux_sign * (_u_secondary[_qp] - _u_primary[_qp]) - eq_pot_qp)/_RTF/2.0));
                }

            }

        }
//        else return _lambda[_qp];
        return residual;
    }

    default:
      return 0;
  }
}


void
GapDisplacementConductanceConstraint::computeConductance(const ADReal & gap)
{
 if (_k_function)
 {
     auto p = _q_point[_qp];
//     const Point p;
//     if (gap > std::numeric_limits<Real>::epsilon()){
//         _k = 0.0;
//         return;
//     }
     _k = _k_function->value(MetaPhysicL::raw_value(gap), p);
//     if (gap > 1e-4) _k = 0.0;
 }
}
