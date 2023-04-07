//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

/**
 * Material objects inherit from Material and override computeQpProperties.
 *
 * Their job is to declare properties for use by other objects in the
 * calculation such as Kernels and BoundaryConditions.
 */
class paperHongli : public Material
{
public:
  static InputParameters validParams();

  paperHongli(const InputParameters & parameters);

protected:
  /// Necessary override. This is where the values of the properties are computed.
  virtual void computeQpProperties() override;

  const ADVariableValue & _cLi;
  const Real & _c1; 
  const Real & _c2;

  /// The diffusivity um^2/s
  ADMaterialProperty<Real> & _diffusivity;
  /// The conductivity mS/cm
  ADMaterialProperty<Real> & _ionic_conductivity;

};
