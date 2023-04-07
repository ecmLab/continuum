/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   EqualNormalFluxConstraint.h
 * Author: srinath
 *
 * Created on January 28, 2021, 1:21 PM
 */

#pragma once

#include "ADMortarConstraint.h"

/**
 * Constrain a specified component of the gradient of a variable to be the same
 * on both sides of an interface.
 */
class EqualNormalFluxConstraint : public ADMortarConstraint
{
public:
  static InputParameters validParams();

  EqualNormalFluxConstraint(const InputParameters & parameters);

protected:
  ADReal computeQpResidual(Moose::MortarType mortar_type) final;

  bool _primary_tensor;
  bool _secondary_tensor;
  
  const ADMaterialProperty<Real> * _primary_mat_prop_real;
  const ADMaterialProperty<Real> * _secondary_mat_prop_real;

  const ADMaterialProperty<RealTensorValue> * _primary_mat_prop_tensor;
  const ADMaterialProperty<RealTensorValue> * _secondary_mat_prop_tensor;
  
};

