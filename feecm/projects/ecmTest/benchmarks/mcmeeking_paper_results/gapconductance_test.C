/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   GapDisplacementConductanceConstraint.C
 * Author: srinath
 *
 * Created on October 3, 2019, 8:16 PM
 */

#include "GapDisplacementConductanceConstraint.h"
#include "Function.h"
#include "SubProblem.h"
#include "MaterialBase.h"

registerADMooseObject("electro_chemo_mechApp", GapDisplacementConductanceConstraint);

InputParameters
GapDisplacementConductanceConstraint::validParams()
{
	InputParameters params = ADMortarConstraint::validParams();
    params.addClassDescription(
        "Computes the residual and Jacobian contributions for the 'Lagrange Multiplier' "
        "implementation of the thermal contact problem. For more information, see the "
        "detailed description here: http://tinyurl.com/gmmhbe9");
    params.addParam<FunctionName>("k_function", "Gap conductance function");
    params.addParam<Real>("k", "Gap conductance");
    params.addCoupledVar("displacements", "Displacement variables");
    params.addParam<NonlinearVariableName>("concentration_primary", "Concentration variable name");
    params.addParam<NonlinearVariableName>("concentration_secondary", "Concentration variable name");
    params.addParam<bool>("include_equilibrium_potential", false, "Equilibrium Potential of Electrode");
    params.addParam<bool>("include_gap", false, "Include gap value in calculation");
    params.addParam<bool>("include_concentration", false, "Include concentration value in calculation");
    MooseEnum computeType("FLUXLINEAR FLUXBUTLERVOMER, EQUILIBRIATE", "FLUXLINEAR");
    params.addParam<MooseEnum>("computeType", computeType, "What do we want this constraint to do");
    params.addParam<bool>("secondary_flux_surface", true, "Default for FLUX approach is secondary surface");
    return params;
}


GapDisplacementConductanceConstraint::GapDisplacementConductanceConstraint(const InputParameters & parameters)
  : ADMortarConstraint(parameters),
    _k_function(isParamValid("k_function") ? & this->getFunction("k_function") : NULL),
    _my_k(getParam<Real>("k")),
    _disp_name(parameters.getVecMooseType("displacements")),
    _n_disp(_disp_name.size()),
    _disp_secondary(_n_disp),
    _disp_primary(_n_disp),
    _include_concentration(getParam<bool>("include_concentration")),
    _concentration_var_primary((isParamValid("concentration_primary") && _include_concentration)
                          ? &this->_subproblem.getStandardVariable(
                                _tid, parameters.getMooseType("concentration_primary"))
                          : nullptr),
    _concentration_var_secondary((isParamValid("concentration_secondary") && _include_concentration)
                          ? &this->_subproblem.getStandardVariable(
                                _tid, parameters.getMooseType("concentration_secondary"))
                          : nullptr),
    _include_equilibrium_potential(getParam<bool>("include_equilibrium_potential")),
    _include_gap(getParam<bool>("include_gap")),
    _compute(getParam<MooseEnum>("computeType").getEnum<Compute>()),
    _secondary_flux_surface(getParam<bool>("secondary_flux_surface"))

{
    if (_k_function == NULL && !parameters.isParamSetByUser("k"))
        mooseError("Either k or k_function must be given");
    if (_concentration_var_primary)
    {
        _concentration_primary = &_concentration_var_secondary->adSlnNeighbor();
    }
      if (_concentration_var_secondary)
    {
        _concentration_secondary = &_concentration_var_secondary->adSln();
    }
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
        _primary_equilibrium_potential = NULL;
        _secondary_equilibrium_potential = NULL;
        _secondary_equilibrium_potential = &getADMaterialProperty<Real>("equilibrium_potential");
        if (&getADMaterialProperty<Real>("equilibrium_potential"))
        {
            _secondary_equilibrium_potential = &getADMaterialProperty<Real>("equilibrium_potential");
        }
        if (_primary_equilibrium_potential && &getNeighborADMaterialProperty<Real>("equilibrium_potential"))
        {
            _primary_equilibrium_potential = &getNeighborADMaterialProperty<Real>("equilibrium_potential");
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
      return _lambda[_qp] * _test_secondary[_i][_qp];
    case Moose::MortarType::Lower:
    {
        bool conc_primary = (_concentration_primary && _include_concentration && ((*_concentration_primary)[_qp] < 0));
        bool conc_secondary (_concentration_secondary && _include_concentration && ((*_concentration_secondary)[_qp] < 0));
        /// Check which concentration variable is available
        bool check_conc = (_secondary_flux_surface ? conc_secondary : conc_primary);
        auto residual = _lambda[_qp] * _test[_i][_qp];
        _k = _my_k;
        if (_has_primary)
        {
            if (_include_gap)
            {
                // we are creating a dual version of phys points primary and secondary here...
                DualRealVectorValue dual_phys_points_primary;
                DualRealVectorValue dual_phys_points_secondary;
                DualRealVectorValue dual_normals;
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
            if  (_compute == Compute::FluxLinear || _compute == Compute::FluxButlerVolmer)
            {
                auto equilibrium_potential = _secondary_equilibrium_potential;
                if (!_secondary_flux_surface ) equilibrium_potential = _primary_equilibrium_potential;
                if (_include_equilibrium_potential && equilibrium_potential)
                {
                    if (check_conc)
                        return residual;
                    else
                    {
                        if (_compute == Compute::FluxLinear)
                        residual = _test[_i][_qp] *
                                (_lambda[_qp] + _k * (_u_secondary[_qp] - _u_primary[_qp] - (*equilibrium_potential)[_qp]));
                        else if (_compute == Compute::FluxButlerVolmer)
                        {
                            /// TBD add exponential terms
                            residual = _test[_i][_qp] *
                                (_lambda[_qp] + _k * (_u_secondary[_qp] - _u_primary[_qp] - (*equilibrium_potential)[_qp]));
                        }
                    }
                }
                else
                {
                    if (check_conc)
                        return residual;
                    else
                    {
                        if (_compute == Compute::FluxLinear)
                            residual = _test[_i][_qp] *
                                (_lambda[_qp] + _k * (_u_secondary[_qp] - _u_primary[_qp]));
                        else if (_compute == Compute::FluxButlerVolmer) {
                            residual = _test[_i][_qp] *
                                (_lambda[_qp] + _k * (_u_secondary[_qp] - _u_primary[_qp]));
                        }
                    }
                }
            }
            else if (_compute == Compute::Equilibriate)
            {
                /// Here we need both concentration variables available
                /// Return 0 flux if either concentration is negative
                if (conc_primary || conc_secondary)
                    return residual;
                /// Return 0 flux if euqilibrium potential is not present on both sides
                if (!_primary_equilibrium_potential || !_secondary_equilibrium_potential)
                    return residual;
                /// check for direction
                residual = _test[_i][_qp] * (_lambda[_qp]
                            + _k * ((*_primary_equilibrium_potential)[_qp] -
                                    (*_secondary_equilibrium_potential)[_qp]));
            }
        }
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
