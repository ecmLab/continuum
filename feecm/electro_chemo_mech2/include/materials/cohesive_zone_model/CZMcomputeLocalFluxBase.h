/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   computeLocalFluxBase.h
 * Author: srinath
 *
 * Created on January 25, 2022, 1:57 PM
 */

#pragma once

#include "InterfaceMaterial.h"
/**
 * Base class used to implement traction separetion laws. All traction separetion laws
 * shall compute the interface traction using the interface coordinate system, and
 * traction derivtives w.r.t. to the interface displacement jump. Interface traction and related
 * derivatives should be implemented overriding the computeInterfaceTractionAndDerivatives method.
 * The interface coordinate system assumes the three component of the traction and
 * disaplcement jump being ordered as [N,S1,S2], where N is the normal component and S1, S2 two
 * orthogonal tangential components. The model also assumes isotropic behavior in the tangential
 * directions.
 */
class CZMcomputeLocalFluxBase : public InterfaceMaterial
{
public:
  static InputParameters validParams();
  CZMcomputeLocalFluxBase(const InputParameters & parameters);

protected:
  void initQpStatefulProperties() override;
  void computeQpProperties() override;

  /// Compute the local flux and derivatives. 
  ///This method should fill the _interface_flux and _dinterface_flux_djump and 
  /// _dinterface_flux_dvariablejump varaibles
  
  virtual void computeInterfaceFluxAndDerivatives() = 0;

  /// Base name of the material system
  const std::string _base_name;
  
  const bool _include_gap;

  /// the value of the flux in local coordinates
  MaterialProperty<Real> & _interface_flux;
  
  /// the flux derivatives wrt the variable in local coordinates
  MaterialProperty<Real> & _dinterface_flux_dvariablejump;

  /// the flux derivatives wrt the displacement jump in local coordinates
  MaterialProperty<RealVectorValue> & _dinterface_flux_djump;

  /// The variable jump in local coordaintes
  const MaterialProperty<Real> & _interface_variable_jump;
  
  const MaterialProperty<RealVectorValue> * _interface_displacement_jump;
  
};
