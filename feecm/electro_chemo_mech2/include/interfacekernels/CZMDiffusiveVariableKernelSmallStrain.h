/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CZMDiffusiveVariableKernelSmallStrain.h
 * Author: srinath
 *
 * Created on January 25, 2022, 8:27 AM
 */

#pragma once
#include "CZMDiffusiveVariableKernelBase.h"

/// DG cohesive zone model kernel for the small strain formulation for diffusive 
/// type variables. 
/// This kernel assumes that the flux is dependent only on the jump of the diffusive
/// variable across the interface. 

class CZMDiffusiveVariableKernelSmallStrain : public CZMDiffusiveVariableKernelBase
{
public:
  static InputParameters validParams();
  CZMDiffusiveVariableKernelSmallStrain(const InputParameters & parameters);

protected:
  Real computeDResidualDVariable(const Moose::DGJacobianType & type) const override;
    
  Real computeDResidualDDisplacement(const unsigned int & component_j,
                                     const Moose::DGJacobianType & type) const override;
  
};


