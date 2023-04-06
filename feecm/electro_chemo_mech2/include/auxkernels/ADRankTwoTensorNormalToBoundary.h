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


class ADRankTwoTensorNormalToBoundary;

// template <>
// InputParameters validParams<ADRankTwoTensorNormalToBoundary>();

/**
 * AD version Auxiliary kernel responsible for computing the components of the flux vector
 * in diffusion problems in a direction normal to the boundary
 */
class ADRankTwoTensorNormalToBoundary : public AuxKernel
{
public:
  static InputParameters validParams();
  ADRankTwoTensorNormalToBoundary(const InputParameters & parameters);

protected:
  virtual Real computeValue();

   /// the normals along the slave face
  const MooseArray<Point> & _normals;
  

  /// Holds the diffusivity from the material system
  const ADMaterialProperty<RankTwoTensor> & _tensor;
  
//  //// Boundary name 
//  BoundaryName _bnd_name;
};

