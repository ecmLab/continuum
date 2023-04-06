/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CZMComputeVariableJumpBase.h
 * Author: srinath
 *
 * Created on January 25, 2022, 1:06 PM
 */

#pragma once

#include "InterfaceMaterial.h"
/**
 * This interface material class computes the jump in the local variable 
 */
class CZMComputeVariableJumpBase : public InterfaceMaterial
{
public:
  static InputParameters validParams();
  CZMComputeVariableJumpBase(const InputParameters & parameters);

protected:
  void computeQpProperties() override;
//  void initQpStatefulProperties() override;

  /// method used to compute the disaplcement jump in interface coordinates according to a
  ///  specific kinematic formulation
  virtual void computeLocalVariableJump() = 0;
  
  /// Base name of the material system
  const std::string _base_name;

  const bool _include_gap;

  /// the rotation matrix transforming from the interface to the global coordinate systems
//  const MaterialProperty<RankTwoTensor> & _czm_total_rotation;
  
  const MaterialProperty<RankTwoTensor> * _czm_total_rotation;
  MaterialProperty<RankTwoTensor> & _total_rotation;
  
  MaterialProperty<Real> & _variable_jump_global;
  MaterialProperty<Real> & _interface_variable_jump;
  
    /// Pointer to the variable whose jump is calculated
  const VariableValue * _variable;
  const VariableValue * _variable_neighbor;
  
};


