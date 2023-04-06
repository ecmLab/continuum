/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   SalehanIrani3DCTractionViscosity.C
 * Author: srinath
 * 
 * Created on January 27, 2022, 6:29 AM
 */

#include "SalehaniIrani3DCTractionViscosity.h"

registerMooseObject("electro_chemo_mechApp", SalehaniIrani3DCTractionViscosity);
    
InputParameters
SalehaniIrani3DCTractionViscosity::validParams()
{
    InputParameters params = SalehaniIrani3DCTraction::validParams();
    params.addClassDescription("3D Coupled (3DC) cohesive law of Salehani and Irani with viscous stabilization");
    params.addRequiredParam<Real>(
      "normal_viscosity",
      "The value of viscosity for normal component");
    params.addRequiredParam<Real>(
      "tangential_viscosity",
      "The value of viscosity for tangential components");
    return params;
}

SalehaniIrani3DCTractionViscosity::SalehaniIrani3DCTractionViscosity(const InputParameters & parameters)
  : SalehaniIrani3DCTraction(parameters),
    _viscosity({getParam<Real>("normal_viscosity"),
               getParam<Real>("tangential_viscosity"),
               getParam<Real>("tangential_viscosity")}), 
    _interface_displacement_jump_inc(
        declarePropertyByName<RealVectorValue>(_base_name + "interface_displacement_jump_inc")),
    _interface_displacement_jump_old(
        getMaterialPropertyOldByName<RealVectorValue>(_base_name + "interface_displacement_jump"))

{
}

void
SalehaniIrani3DCTractionViscosity::computeInterfaceTractionAndDerivatives()
{
   _interface_displacement_jump_inc[_qp] =
        _interface_displacement_jump[_qp] - _interface_displacement_jump_old[_qp];
  _interface_traction[_qp] = computeTraction();
  _dinterface_traction_djump[_qp] = computeTractionDerivatives();
}

RealVectorValue
SalehaniIrani3DCTractionViscosity::computeTraction()
{
      // The convention for ordering the traction is N, T, S, where N is the normal direction, and T and
  // S are two arbitrary tangential directions.
  RealVectorValue traction_local;
    traction_local = SalehaniIrani3DCTraction::computeTraction();
    for (unsigned int i = 0; i<3; i++)
    {
        traction_local(i) *= (1.0 + _viscosity(i) / _max_allowable_traction(i) * 
                            _interface_displacement_jump_inc[_qp](i)/_delta_u0(i)/_dt);
    }
    return traction_local;
}

RankTwoTensor
SalehaniIrani3DCTractionViscosity::computeTractionDerivatives()
{
    RankTwoTensor traction_jump_derivatives_local;
    traction_jump_derivatives_local = SalehaniIrani3DCTraction::computeTractionDerivatives();
    for (unsigned int i = 0; i <3; i++) 
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            traction_jump_derivatives_local(i,j) *= (1.0 + _viscosity(i)/ _max_allowable_traction(i) * 
                    _interface_displacement_jump_inc[_qp](i)/_delta_u0(i)/_dt);
            if (i == j)
            {
                traction_jump_derivatives_local(i,i)  += _interface_traction[_qp](i) * 
                _viscosity(i)/_delta_u0(i)/_dt/_max_allowable_traction(i);
            }
            
        }
    }
    return traction_jump_derivatives_local;
}