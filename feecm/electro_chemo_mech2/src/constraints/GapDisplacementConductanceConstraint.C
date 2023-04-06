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
    return params;
}


GapDisplacementConductanceConstraint::GapDisplacementConductanceConstraint(const InputParameters & parameters)
  : ADMortarConstraint(parameters),
        _k_function(isParamValid("k_function") ? & this->getFunction("k_function") : NULL),
        _my_k(getParam<Real>("k")), 

        _disp_name(parameters.getVecMooseType("displacements")),
     _n_disp(_disp_name.size()),
     _disp_secondary(_n_disp),
     _disp_primary(_n_disp)
{
//    if (_k_function != NULL && parameters.isParamSetByUser("k"))
//        mooseError("Both k and k_function cannot be given at the same time");
    if (_k_function == NULL && !parameters.isParamSetByUser("k"))
        mooseError("Either k or k_function must be given");

    for (unsigned int i = 0; i < _n_disp; ++i)
    {
        auto & disp_var = this->subProblem().getStandardVariable(_tid, _disp_name[i]);
        _disp_secondary[i] = &disp_var.adSln();
        _disp_primary[i] = &disp_var.adSlnNeighbor();
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
//      auto l = (_phys_points_primary[_qp] - _phys_points_secondary[_qp]).norm();
      // Normal distance across the gap
//      auto l = (_phys_points_primary[_qp] - _phys_points_secondary[_qp]) * _normals[_qp];
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
  //        dual_normals(i).derivatives() = _normals[_qp](i).derivatives();
        }

  //      auto l = (dual_phys_points_primary - dual_phys_points_secondary).norm();
  //      auto l = (_phys_points_primary[_qp] - _phys_points_secondary[_qp]) * _normals[_qp];
        auto l = (dual_phys_points_primary - dual_phys_points_secondary) * dual_normals;
      if (_k_function)
        computeConductance(l);
      else
          _k = _my_k;
      return (_k * (_u_primary[_qp] - _u_secondary[_qp]) - _lambda[_qp]) * _test[_i][_qp];
    } 
      
    default:
      return 0;
  }
}

//template <>
//DualReal
//GapDisplacementConductanceConstraint<JACOBIAN>::computeQpResidual(Moose::MortarType mortar_type)
//{
//  switch (mortar_type)
//  {
//    case Moose::MortarType::Primary:
//      return _lambda[_qp] * _test_primary[_i][_qp];
//
//    case Moose::MortarType::Secondary:
//      return -_lambda[_qp] * _test_secondary[_i][_qp];
//
//    case Moose::MortarType::Lower:
//    {
//      // we are creating a dual version of phys points primary and secondary here...
//      DualRealVectorValue dual_phys_points_primary;
//      DualRealVectorValue dual_phys_points_secondary;
//      DualRealVectorValue dual_normals;
//      for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
//      {
//        dual_phys_points_primary(i) = _phys_points_primary[_qp](i);
//        dual_phys_points_secondary(i) = _phys_points_secondary[_qp](i);
//        dual_normals(i) = _normals[_qp](i);
//      }
//
//      // ...which uses the derivative vector of the primary and secondary displacements as
//      // an approximation of the true phys points derivatives
//      for (unsigned int i = 0; i < _n_disp; ++i)
//      {
//        dual_phys_points_primary(i).derivatives() = (*_disp_primary[i])[_qp].derivatives();
//        dual_phys_points_secondary(i).derivatives() = (*_disp_secondary[i])[_qp].derivatives();
////        dual_normals(i).derivatives() = _normals[_qp](i).derivatives();
//      }
//
////      auto l = (dual_phys_points_primary - dual_phys_points_secondary).norm();
////      auto l = (_phys_points_primary[_qp] - _phys_points_secondary[_qp]) * _normals[_qp];
//      auto l = (dual_phys_points_primary - dual_phys_points_secondary) * dual_normals;
//      if (_k_function)
//        computeConductance(l);
//      else
//          _k = _my_k;
//        return (_k * (_u_primary[_qp] - _u_secondary[_qp]) - _lambda[_qp]) * _test[_i][_qp];
//    }
//
//    default:
//      return 0;
//  }
//}


void 
GapDisplacementConductanceConstraint::computeConductance(const ADReal & penetration)
{
 if (_k_function)
 {
     const Point p;
     _k = _k_function->value(MetaPhysicL::raw_value(penetration), p);
//     if (penetration > 1e-8) _k = 0.0;
 }
}


//template<>
//void
//GapDisplacementConductanceConstraint<JACOBIAN>::computeConductance(const DualReal & penetration)
//{
//    if (_k_function)
//    {
//        const Point p;
//        _k = _k_function->value(MetaPhysicL::raw_value(penetration), p);
//        _k.derivatives() = penetration.derivatives()*
//                _k_function->timeDerivative(MetaPhysicL::raw_value(penetration),p);
//    }
//}
