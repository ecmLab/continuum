/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   ADDiffusionFluxAux.h
 * Author: srinath
 *
 * Created on October 18, 2019, 8:50 AM
 */


#pragma once

#include "AuxKernel.h"

class ADDiffusionFluxAux;

// template <>
// InputParameters validParams<ADDiffusionFluxAux>();

/**
 * Auxiliary kernel responsible for computing the components of the flux vector
 * in diffusion problems
 */
class ADDiffusionFluxAux : public AuxKernel
{
public:
  static InputParameters validParams();

  ADDiffusionFluxAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();
  /// Will hold 0, 1, or 2 corresponding to x, y, or z.
  int _component;

  /// Holds the solution gradient at the current quadrature points
  const VariableGradient & _grad_u;

  /// Holds the diffusivity from the material system
  const ADMaterialProperty<Real> & _diffusion_coef;
};

