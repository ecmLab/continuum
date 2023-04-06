/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   czm_variable_base.h
 * Author: srinath
 *
 * Created on January 24, 2022, 11:06 AM
 */

#pragma once

#include "InterfaceKernel.h"
#include "JvarMapInterface.h"

/// Base class for implementing DG cohesive zone models (CZM) for 1D,2D, and 3D
/// separation laws for any diffusive type variable other than the displacements, 
/// eg. voltage, concentration, temperature etc. 
/// This kernel operates on the variable and also computes the off diagonal
/// jacobian contributions to the displacements. This necessitates
/// the existence of a displacement based cohesive zone model and the relevant 
/// material properties. The single variable will only have the normal component
/// of the flux that is affected by the cohesive zone and only the normal 
//// component of the displacement

class CZMDiffusiveVariableKernelBase : public JvarMapKernelInterface<InterfaceKernel>
{
public:
  static InputParameters validParams();
  CZMDiffusiveVariableKernelBase(const InputParameters & parameters);
 protected:
  Real computeQpResidual(Moose::DGResidualType type) override;
  Real computeQpJacobian(Moose::DGJacobianType type) override;
  Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;
  
  // method computing derivative of residual w.r.t variable
  virtual Real computeDResidualDVariable(const Moose::DGJacobianType & type) const = 0;
  
  /// method computing the derivative of residual[_component] w.r.t displacement[component_j]
  virtual Real computeDResidualDDisplacement(const unsigned int & component_j,
                                             const Moose::DGJacobianType & type) const = 0;

  /// Base name of the material system that this kernel applies to
  const std::string _base_name;

  const bool _include_gap;
  
  /// the displacement component this kernel is operating on (0=x, 1=y, 2 =z)
  const unsigned int _component;

  /// number of displacement components
  const unsigned int _ndisp;

  /// Coupled displacement component variable IDs
  ///@{
  std::vector<unsigned int> _disp_var;
  std::vector<unsigned int> _disp_neighbor_var;
  ///@}

  // pointer to displacement variables
  std::vector<MooseVariable *> _disp_vars;
  

//  /// Coupled variable ID
//  const unsigned int _var_id;
//  
//  /// Pointer to the coupled variable
//  MooseVariable * _var;
  
  // values of the flux and flux derivatives used
  ///@{
  const MaterialProperty<Real> & _flux_global;
  const MaterialProperty<Real> & _dflux_dvariablejump_global;
  const MaterialProperty<RealVectorValue> & _dflux_djump_global;
  ///@}
  
};


