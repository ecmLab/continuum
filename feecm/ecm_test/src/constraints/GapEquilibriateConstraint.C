/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   GapEquilibriateConstraint.C
 * Author: srinath
 * 
 * Created on January 19, 2021, 3:50 PM
 */

#include "GapEquilibriateConstraint.h"

#include "Function.h"
#include "SubProblem.h"
#include "MaterialBase.h"

registerADMooseObject("electro_chemo_mechApp", GapEquilibriateConstraint);

InputParameters
GapEquilibriateConstraint::validParams()
{
	InputParameters params = ADMortarConstraint::validParams();
    params.addClassDescription(
        "Computes the residual and Jacobian contributions for the 'Lagrange Multiplier' "
        "implementation of the thermal contact problem. For more information, see the "
        "detailed description here: http://tinyurl.com/gmmhbe9");
    params.addParam<FunctionName>("k_function", "Function of the relevant gap problem");
    params.addParam<Real>("k", "Gap coefficient");
    params.addParam<Real>("faraday", 96485.3329, "Faraday's Constant");
    params.addParam<Real>("R", 8.3145, "Universal Gas Constant");
    params.addParam<Real>("temperature", 298, "Value of temperature to use");
    params.addCoupledVar("displacements", "Displacement variables");
    params.addParam<bool>("include_gap", false, "Include gap value in calculation");
      params.addParam<MaterialPropertyName>(
      "primary_mat_prop",
      "equilibrium_potential",
      "The material property name providing the quantity to equilibriate on the primary side");
    params.addParam<MaterialPropertyName>(
      "secondary_mat_prop",
      "equilibrium_potential",
      "The material property name providing the quantity to equilibriate on the secondary side");
    MooseEnum oneSidedTransfer("PRIMARY->SECONDARY SECONDARY->PRIMARY NONE", "NONE");
    params.addParam<MooseEnum>("one_sided", oneSidedTransfer, "Is the equilibriation one sided");
    MooseEnum preFactor("RT/F RT NONE", "RT/F");
    params.addParam<MooseEnum>("prefactor", preFactor, "Pre-factor to divide material property");
    return params;
}


GapEquilibriateConstraint::GapEquilibriateConstraint(const InputParameters & parameters)
  : ADMortarConstraint(parameters),
    _k_function(isParamValid("k_function") ? & this->getFunction("k_function") : NULL),
    _my_k(getParam<Real>("k")), 
    _gas_constant(getParam<Real>("R")),
    _temp(getParam<Real>("temperature")),
    _faraday(getParam<Real>("faraday")),
    _disp_name(parameters.getVecMooseType("displacements")),
    _n_disp(_disp_name.size()),
    _disp_secondary(_n_disp),
    _disp_primary(_n_disp), 
    _include_gap(getParam<bool>("include_gap")),
    _one_sided(getParam<MooseEnum>("one_sided").getEnum<OneSided>()),
    _pref(getParam<MooseEnum>("prefactor").getEnum<preFactor>())
{
    if (_k_function == NULL && !parameters.isParamSetByUser("k"))
        mooseError("Either k or k_function must be given");
    if (_include_gap)
    {
        for (unsigned int i = 0; i < _n_disp; ++i)
        {
            auto disp_var = &this->_subproblem.getStandardVariable(_tid, _disp_name[i]);
            _disp_secondary[i] = &disp_var->adSln();
            _disp_primary[i] = &disp_var->adSlnNeighbor();
        }
    }
    _primary_mat_prop = (&getNeighborADMaterialProperty<Real>("primary_mat_prop") ?
                                        &getNeighborADMaterialProperty<Real>("primary_mat_prop") : nullptr);
    _secondary_mat_prop = (&getADMaterialProperty<Real>("secondary_mat_prop") ?
                                        &getADMaterialProperty<Real>("secondary_mat_prop") : nullptr);
    if (_pref == preFactor::RTF)
        _prefactor = _gas_constant * _temp / _faraday;
    else if (_pref == preFactor::RT)
        _prefactor = _gas_constant * _temp;
    else if (_pref == preFactor::None)
        _prefactor = 1.0;
}

ADReal
GapEquilibriateConstraint::computeQpResidual(Moose::MortarType mortar_type)
{
    
  switch (mortar_type)
  {
    case Moose::MortarType::Primary:
      return _lambda[_qp] * _test_primary[_i][_qp];
    case Moose::MortarType::Secondary:
      return _lambda[_qp] * _test_secondary[_i][_qp];
    case Moose::MortarType::Lower:
    {
        auto residual = _lambda[_qp] * _test[_i][_qp];
        _k = _my_k;
//        if (_has_primary) 
        {
            if (_include_gap)
            {
                // we are creating a dual version of phys points primary and secondary here...
                RealVectorValue dual_phys_points_primary;
                RealVectorValue dual_phys_points_secondary;
                RealVectorValue dual_normals;
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
            
            /// (TBD) check to make sure sign is right when primary and secondary
            /// surfaces are switched .... 
//            if (_u_primary[_qp] < 0 || _u_secondary[_qp] < 0)
//                return residual;
            if (!_primary_mat_prop || !_secondary_mat_prop)
                return residual;
            
            auto diff_mat_prop_qp = _k/_prefactor * ((*_primary_mat_prop)[_qp] - (*_secondary_mat_prop)[_qp]);
//            auto diff_mat_prop_qp = _k/_prefactor * (_u_primary[_qp] - _u_secondary[_qp]);
            if (_one_sided == OneSided::Primary_Secondary)
            {
                if (diff_mat_prop_qp > 0) diff_mat_prop_qp = 0.0;
            }
            if (_one_sided == OneSided::Secondary_Primary)
            {
                if (diff_mat_prop_qp < 0) diff_mat_prop_qp = 0.0;
            }
            
            residual = _test[_i][_qp] * 
                    (_lambda[_qp] + diff_mat_prop_qp);
        }
        return residual;
    } 
      
    default:
      return 0;
  }
}


void 
GapEquilibriateConstraint::computeConductance(const ADReal & gap)
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

