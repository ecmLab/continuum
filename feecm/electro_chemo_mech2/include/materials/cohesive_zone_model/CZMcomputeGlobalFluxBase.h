/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   computeGlobalFlux.h
 * Author: srinath
 *
 * Created on January 25, 2022, 2:09 PM
 */


#pragma once

#include "InterfaceMaterial.h"
/**
 * Base class flux computing the flux used to impose equilibrium and its derivatives  w.r.t.
 * the global displacement jump, starting from the values provided from any CZM constituive material
 * model.
 */
class CZMcomputeGlobalFluxBase : public InterfaceMaterial
{
public:
  static InputParameters validParams();
  CZMcomputeGlobalFluxBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// method computing the equilibrium flux and its derivatives
  virtual void computeEquilibriumFluxAndDerivatives() = 0;

  /// Base name of the material system
  const std::string _base_name;

  const bool _include_gap;
  /// the value of the flux in global and interface coordinates
  ///@{
  MaterialProperty<Real> & _flux_global;
  const MaterialProperty<Real> & _interface_flux;
  ///@}

  /// the flux's derivatives w.r.t. the displacement jump in global and interface coordinates
  ///@{
  MaterialProperty<Real> & _dflux_dvariablejump_global;
  MaterialProperty<RealVectorValue> & _dflux_djump_global;
  const MaterialProperty<Real> & _dinterface_flux_dvariablejump;
  const MaterialProperty<RealVectorValue> & _dinterface_flux_djump;
  ///@}

  /// the rotation matrix trnasforming from interface to global coordinates
  const MaterialProperty<RankTwoTensor> & _total_rotation;
};

