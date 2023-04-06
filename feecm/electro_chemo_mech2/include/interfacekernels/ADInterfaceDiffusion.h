/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   ADInterfaceDiffusion.h
 * Author: srinath
 *
 * Created on March 28, 2022, 5:01 PM
 */

#pragma once

#include "ADInterfaceKernel.h"

/**
 * DG kernel for interfacing diffusion between two variables on adjacent blocks
 */
class ADInterfaceDiffusion : public ADInterfaceKernel
{
public:
  static InputParameters validParams();

 ADInterfaceDiffusion(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual(Moose::DGResidualType type) override;
//  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

  const ADMaterialProperty<Real> & _D;
  const ADMaterialProperty<Real> & _D_neighbor;
};
