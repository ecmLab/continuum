/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   AnsioDiffusionFluxNormaltoBoundaryAux.h
 * Author: srinath
 *
 * Created on October 18, 2019, 8:50 AM
 */

#pragma once

#include "AuxKernel.h"


class AnisoDiffusionFluxNormalToBoundary;

// template <>
// InputParameters validParams<AnisoDiffusionFluxNormalToBoundary>();

/**
 * AD version Auxiliary kernel responsible for computing the components of the flux vector
 * in diffusion problems in a direction normal to the boundary
 */
class AnisoDiffusionFluxNormalToBoundary : public AuxKernel
{
public:
  static InputParameters validParams();
  AnisoDiffusionFluxNormalToBoundary(const InputParameters & parameters);

protected:
  virtual Real computeValue();

   /// the normals along the slave face
  const MooseArray<Point> & _normals;

  /// Holds the solution gradient at the current quadrature points
  const VariableGradient & _grad_u;
  

  /// Holds the diffusivity from the material system
  const ADMaterialProperty<RealTensorValue> & _diffusion_coef;
};

