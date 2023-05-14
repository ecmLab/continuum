/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CZMComputeVariableJumpSmallStrain.h
 * Author: srinath
 *
 * Created on January 25, 2022, 1:49 PM
 */

#pragma once

#include "CZMComputeVariableJumpBase.h"
/**
 * Compute the interface displacement jump across a cohesive zone under the small strain
 * assumption
 */
class CZMComputeVariableJumpSmallStrain : public CZMComputeVariableJumpBase
{
public:
  static InputParameters validParams();
  CZMComputeVariableJumpSmallStrain(const InputParameters & parameters);

protected:
  /// compute the total displacement jump in interface coordinates
  void computeLocalVariableJump() override;
};
