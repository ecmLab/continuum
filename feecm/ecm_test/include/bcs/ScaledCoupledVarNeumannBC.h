/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ScaledCoupledVarNeumannBC.h
 * Author: srinath
 *
 * Created on October 18, 2019, 3:01 PM
 */

#pragma once

#include "IntegratedBC.h"

class ScaledCoupledVarNeumannBC;

// template <>
// InputParameters validParams<ScaledCoupledVarNeumannBC>();

/**
 * Implements a Neumann BC where grad(u)=_coupled_var on the boundary.
 * Uses the term produced from integrating the diffusion operator by parts.
 */
class ScaledCoupledVarNeumannBC : public IntegratedBC
{
public:

  static InputParameters validParams();
  ScaledCoupledVarNeumannBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

    /// Scaling factor Neumann BC
  const Real _scale;

  /// Variable providing the value of grad(u) on the boundary.
  const VariableValue & _coupled_var;
  
};

