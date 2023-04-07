/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   PressureLMConductanceConstraint.C
 * Author: srinath
 * 
 * Created on October 3, 2019, 8:16 PM
 */

#include "PressureLMConductanceConstraint.h"
#include "Function.h"

registerADMooseObject("electro_chemo_mechApp", PressureLMConductanceConstraint);

InputParameters
PressureLMConductanceConstraint::validParams()
{
	InputParameters params = ADMortarConstraint::validParams();
    params.addClassDescription(
        "Computes the residual and Jacobian contributions for the 'Lagrange Multiplier' "
        "implementation of the thermal contact problem. For more information, see the "
        "detailed description here: http://tinyurl.com/gmmhbe9");
    params.addParam<Real>("k", "Gap conductance");
    params.addParam<Real>("threshold_pressure", 1e-15, "Value of pressure threshold");
    params.addRequiredParam<NonlinearVariableName>(
      "contact_pressure",
      "The normal contact pressure; oftentimes this may be a separate lagrange multiplier "
      "variable");
    return params;
}


PressureLMConductanceConstraint::PressureLMConductanceConstraint(const InputParameters & parameters)
  : ADMortarConstraint(parameters),
    _my_k(getParam<Real>("k")), 
    _threshold_pressure(getParam<Real>("threshold_pressure")),
    _contact_pressure_name(parameters.getMooseType("contact_pressure"))
{
    auto pressure_var = &this->_subproblem.getStandardVariable(_tid, _contact_pressure_name);
    _contact_pressure = &pressure_var->adSlnLower();
}

ADReal
PressureLMConductanceConstraint::computeQpResidual(Moose::MortarType mortar_type)
{
  switch (mortar_type)
  {
    case Moose::MortarType::Primary:
      return _lambda[_qp] * _test_primary[_i][_qp];
    case Moose::MortarType::Secondary:
      return -_lambda[_qp] * _test_secondary[_i][_qp];
    case Moose::MortarType::Lower:
    {
//      if (_has_primary)
      {
          auto pr = (*_contact_pressure)[_qp];
          if (pr > _threshold_pressure)
          {
              _k = _my_k;
          } else
              _k = 0.0;
        return (_k * (_u_primary[_qp] - _u_secondary[_qp]) - _lambda[_qp]) * _test[_i][_qp];
      }
//      else
//          return _lambda[_qp] * _test[_i][_qp];
//      
    } 
      
    default:
      return 0;
  }
}
