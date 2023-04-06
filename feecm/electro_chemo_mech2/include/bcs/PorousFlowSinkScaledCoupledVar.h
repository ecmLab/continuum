/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   PorousFlowSinkScaledCoupledVar.h
 * Author: srinath
 *
 * Created on October 19, 2021, 6:36 PM
 */
#pragma once
#include "PorousFlowSink.h"

class PorousFlowSinkScaledCoupledVar: public PorousFlowSink
{
public:
  static InputParameters validParams();
  PorousFlowSinkScaledCoupledVar(const InputParameters & parameters);
protected:
  virtual Real computeQpResidual() override; 

   /// Variable providing the value of grad(u) on the boundary.
  const VariableValue & _coupled_var;
    /// Scale factor
      /// Scaling factor Neumann BC
  const Real _scale;

};



