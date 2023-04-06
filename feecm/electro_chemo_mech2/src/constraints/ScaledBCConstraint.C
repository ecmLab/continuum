/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ScaledBCConstraint.C
 * Author: srinath
 * 
 * Created on September 18, 2020, 6:38 AM
 */

#include "ScaledBCConstraint.h"

registerADMooseObject("electro_chemo_mechApp", ScaledBCConstraint);

InputParameters
ScaledBCConstraint::validParams()
{
	InputParameters params = ADMortarConstraint::validParams();
    params.addClassDescription("Scales an existing contact problem variable"
                               "to apply on another variable");
    params.addRequiredParam<Real>("scale", "Scaling of flux");
    params.addParam<bool>("primary", true, "Apply bc on primary surface");
    params.set<bool>("compute_lm_residuals") = false;
    
    return params;
}

ScaledBCConstraint::ScaledBCConstraint(const InputParameters & parameters)
  : ADMortarConstraint(parameters), _scale(getParam<Real>("scale")),
        _primary(getParam<bool>("primary"))
{   
}

ADReal
ScaledBCConstraint::computeQpResidual(Moose::MortarType mortar_type)
{
  switch (mortar_type)
  {
    case Moose::MortarType::Primary:
    {
        if (_primary)
            return -_lambda[_qp] * _test_primary[_i][_qp] * _scale;
    }
    case Moose::MortarType::Secondary:
    {
        if (!_primary)
            return -_lambda[_qp] * _test_secondary[_i][_qp] * _scale;
    }
    default:
      return 0;
  }
}
