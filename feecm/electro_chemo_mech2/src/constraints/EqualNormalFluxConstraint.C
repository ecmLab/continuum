/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   EqualNormalFluxConstraint.C
 * Author: srinath
 * 
 * Created on January 28, 2021, 1:21 PM
 */

#include "EqualNormalFluxConstraint.h"

#include "SubProblem.h"
#include "FEProblem.h"

registerMooseObject("electro_chemo_mechApp", EqualNormalFluxConstraint);

InputParameters
EqualNormalFluxConstraint::validParams()
{
  InputParameters params = ADMortarConstraint::validParams();
  params.addClassDescription(
      "EqualNormalFluxConstraint enforces continuity of a gradient component between secondary and "
      "primary sides of a mortar interface using lagrange multipliers");
  params.addParam<MaterialPropertyName>(
      "primary_mat_prop",
      "diffusivity",
      "The material property name providing the quantity to equilibrate on the primary side");
  params.addParam<MaterialPropertyName>(
      "secondary_mat_prop",
      "diffusivity",
      "The material property name providing the quantity to equilibrate on the secondary side");
  params.addParam<bool>("primary_tensor", false, "Is the material property tensor_valued");
  params.addParam<bool>("secondary_tensor", false, "Is the material property tensor_valued");
  
  return params;
}

EqualNormalFluxConstraint::EqualNormalFluxConstraint(const InputParameters & parameters)
  : ADMortarConstraint(parameters),
  _primary_tensor(getParam<bool>("primary_tensor")),
  _secondary_tensor(getParam<bool>("secondary_tensor")),
  _primary_mat_prop_real(!_primary_tensor ? &getNeighborADMaterialProperty<Real>("primary_mat_prop"):nullptr),
  _secondary_mat_prop_real(!_secondary_tensor ? &getADMaterialProperty<Real>("secondary_mat_prop"):nullptr),
  _primary_mat_prop_tensor(_primary_tensor ? &getNeighborADMaterialProperty<RealTensorValue>("primary_mat_prop"):nullptr),
  _secondary_mat_prop_tensor(_secondary_tensor ? &getADMaterialProperty<RealTensorValue>("secondary_mat_prop"):nullptr)

{
}

ADReal
EqualNormalFluxConstraint::computeQpResidual(Moose::MortarType mortar_type)
{
       
  switch (mortar_type)
  {
    case Moose::MortarType::Secondary:
      return _lambda[_qp] * _test_secondary[_i][_qp]; 
    case Moose::MortarType::Primary:
      return _lambda[_qp] * _test_primary[_i][_qp];
    case Moose::MortarType::Lower:
    {
        ADRealVectorValue gradup; 
        ADRealVectorValue gradus;
        auto residual = _lambda[_qp] * _test[_i][_qp];
//        if (_has_primary)
        {
            if (!_primary_tensor)
                gradup = (*_primary_mat_prop_real)[_qp] * _grad_u_primary[_qp];
            else
                gradup = (*_primary_mat_prop_tensor)[_qp] * _grad_u_primary[_qp];
            if (!_secondary_tensor)
                gradus = (*_secondary_mat_prop_real)[_qp] * _grad_u_secondary[_qp];
            else
                gradus = (*_secondary_mat_prop_tensor)[_qp] * _grad_u_secondary[_qp];
                
            residual = ((gradup - gradus) * _normals[_qp]  + _lambda[_qp]) * _test[_i][_qp];
        }
                    
//        return residual;
    }
    default:
      return 0;
  }
}
