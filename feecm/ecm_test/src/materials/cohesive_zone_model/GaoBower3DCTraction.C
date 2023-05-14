/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   GaoBower3DCTraction.C
 * Author: srinath
 * 
 * Created on January 27, 2022, 12:09 PM
 */

#include "GaoBower3DCTraction.h"

registerMooseObject("electro_chemo_mechApp", GaoBower3DCTraction);
    
InputParameters
GaoBower3DCTraction::validParams()
{
    InputParameters params = SalehaniIrani3DCTraction::validParams();
    params.addClassDescription("3D Coupled (3DC) cohesive law of Salehani and Irani with viscous stabilization");
    params.addRequiredParam<Real>(
      "normal_viscosity",
      "The value of viscosity for normal component");
    params.addRequiredParam<Real>(
      "tangential_viscosity",
      "The value of viscosity for tangential components");
    params.addParam<Real>("r", 0, "r parameter");
    params.addParam<Real>("q", 1, "q parameter");
    return params;
}

GaoBower3DCTraction::GaoBower3DCTraction(const InputParameters & parameters)
  : SalehaniIrani3DCTraction(parameters),
    _viscosity({getParam<Real>("normal_viscosity"),
               getParam<Real>("tangential_viscosity"),
               getParam<Real>("tangential_viscosity")}), 
    _r(getParam<Real>("r")),
    _q(getParam<Real>("q")),
    _interface_displacement_jump_inc(
        declarePropertyByName<RealVectorValue>(_base_name + "interface_displacement_jump_inc")),
    _interface_displacement_jump_old(
        getMaterialPropertyOldByName<RealVectorValue>(_base_name + "interface_displacement_jump"))

{
}

void
GaoBower3DCTraction::computeInterfaceTractionAndDerivatives()
{
   _interface_displacement_jump_inc[_qp] =
        _interface_displacement_jump[_qp] - _interface_displacement_jump_old[_qp];
  _interface_traction[_qp] = computeTraction();
  _dinterface_traction_djump[_qp] = computeTractionDerivatives();
}

RealVectorValue
GaoBower3DCTraction::computeTraction()
{
      // The convention for ordering the traction is N, T, S, where N is the normal direction, and T and
  // S are two arbitrary tangential directions.
  RealVectorValue traction_local;
  RealVectorValue normalized_jump, normalized_jump_old, normalized_jump_velocity; 
  for (unsigned int i = 0; i<3; i++)
  {
      normalized_jump(i) = _interface_displacement_jump[_qp](i)/_delta_u0(i);
      normalized_jump_old(i) = _interface_displacement_jump_old[_qp](i)/_delta_u0(i);
      normalized_jump_velocity(i) = (normalized_jump(i)-normalized_jump_old(i))/_dt;
  }
  
  Real exp_x, x, aa;
  x = 0;
  aa = 0;
  Real temp1 = std::exp(1 - normalized_jump(0));
  Real temp2 = _delta_u0(0)/_delta_u0(1);
  for (unsigned int i=1; i<3; i++)
  {
    aa = _interface_displacement_jump[_qp](i) / _delta_u0(i);
    aa *= aa;
    x += aa;   
  }
  exp_x = std::exp(-x);
  
  
          
  /// Xu-Needleman 
  Real bb = 1;
  for (unsigned int i = 0; i<3; i++)
  {
      if (i > 0)
      {
          bb = 2;
          traction_local(i) = bb * temp2 * normalized_jump(i) * 
                  (_q + ((_r-_q)/(_r-1)) * normalized_jump(0)) * (exp_x * temp1);
          
      } else
      {
          bb = 1;
          traction_local(i) =bb * temp1 * (normalized_jump(i) 
                  * exp_x + ((1 - _q)/(_r - 1)) * (1-exp_x) * (_r - normalized_jump(i)));
      }
      traction_local(i) *= _max_allowable_traction(0);
              
  }
    /// Viscosity contribution
    for (unsigned int i = 0; i<3; i++)
    {
        traction_local(i) *= (1.0 + _viscosity(i) / _max_allowable_traction(0) * 
                            _interface_displacement_jump_inc[_qp](i)/_delta_u0(i)/_dt);
    }
    return traction_local;
}

RankTwoTensor
GaoBower3DCTraction::computeTractionDerivatives()
{
    RankTwoTensor traction_jump_derivatives_local;
    unsigned int i,j;
    Real aa, bb, exp_x, x;
    
    RealVectorValue normalized_jump, normalized_jump_old, normalized_jump_velocity; 
    for (i = 0; i<3; i++)
    {
        normalized_jump(i) = _interface_displacement_jump[_qp](i)/_delta_u0(i);
        normalized_jump_old(i) = _interface_displacement_jump_old[_qp](i)/_delta_u0(i);
        normalized_jump_velocity(i) = (normalized_jump(i)-normalized_jump_old(i))/_dt;
    }

    Real temp1 = std::exp(- normalized_jump(0));
    Real exp1 = std::exp(1.0);
    Real temp2 = _delta_u0(0)/_delta_u0(1);
    
    x = 0;
    aa = 0;

    
    for (i=1; i<3; i++)
    {
        aa = _interface_displacement_jump[_qp](i) / _delta_u0(i);
        aa *= aa;
        x += aa;   
    }
    exp_x = std::exp(-x);
    
    for (i = 0; i <3; i++)
    {
        for (j = 0; j<3;j++)
            if (i == 0){
                if (j == 0){
                    traction_jump_derivatives_local(i,j) = 
                            -temp1 * (normalized_jump(j) * exp_x * 
                            ((1-_q)/(_r-1)) * (1 - exp_x)*(_r-normalized_jump(j)));
                    traction_jump_derivatives_local(i,j) += exp1 * temp1 * 
                            (exp_x - ((1-_q)/(_r-1)) * (1 - exp_x));
                } else
                {
                    traction_jump_derivatives_local(i,j) = 2 * exp1 * temp1 * exp_x 
                            * normalized_jump(j);
                    traction_jump_derivatives_local(i,j) *= (- normalized_jump(i) 
                            + (1-_q)/(_r-1) * (_r - normalized_jump(j)));
                            
                }
                traction_jump_derivatives_local(i,j) *= _max_allowable_traction(0) /
                        _delta_u0(j);
            } else
            {
                if (j == 0){
                    traction_jump_derivatives_local(i,j) = 2 * temp2 * normalized_jump(i) * exp_x * exp1 * temp1;
                    traction_jump_derivatives_local(i,j) *= ((_r - _q)/(_r -1) * (1- normalized_jump(j) - _q));                    
                    
//                    traction_jump_derivatives_local(i,j) = 2 * temp2 * normalized_jump(i)
//                            * ((_r - _q)/(_r -1)) * exp1 * temp1 * exp_x;
//                    traction_jump_derivatives_local(i,j) -= 2 * temp2 * normalized_jump(i)
//                            * (_q + ((_r - _q)/(_r -1)) * normalized_jump(0)) * temp1 * exp_x;
                } else
                {
                    traction_jump_derivatives_local(i,j) = 2 * temp2 * exp1 * temp1 * exp_x 
                            * (_q + ((_r - _q)/(_r -1)) * normalized_jump(0));
                    traction_jump_derivatives_local(i,j) *= (1 - 2 * normalized_jump(j)
                            * normalized_jump(j));
                            
//                    traction_jump_derivatives_local(i,j) = 2 * temp2 
//                            * (_q + ((_r - _q)/(_r -1)) * normalized_jump(0)) 
//                            * exp1 * temp1 * exp_x;
//                    
//                    traction_jump_derivatives_local(i,j) -= 2 * temp2 * normalized_jump(i) 
//                            * (_q + ((_r - _q)/(_r -1)) * normalized_jump(0)) 
//                            * exp1 * temp1 * 2*normalized_jump(j);
                }
                traction_jump_derivatives_local(i,j) *= _max_allowable_traction(0) /
                        _delta_u0(j);
                    
            }
            
    }    
    /// Viscosity contribution
    for (i = 0; i <3; i++) 
    {
        for (j = 0; j < 3; j++)
        {
            traction_jump_derivatives_local(i,j) *= (1.0 + _viscosity(i)/ _max_allowable_traction(0) * 
                    normalized_jump_velocity(i));
            if (i == j)
            {
                traction_jump_derivatives_local(i,i)  += _interface_traction[_qp](i) * 
                _viscosity(i)/_delta_u0(i)/_dt/_max_allowable_traction(0);
            }
            
        }
    }
    return traction_jump_derivatives_local;
}