    //* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
/* 
 * File:   ADStressBasedChemicalPotential.h
 * Author: srinath
 * 
 * Created on October 20, 2020, 06:50 AM
 */

#pragma once

#include "ADKernel.h"

class ADStressBasedChemicalPotential: public ADKernel
{
public:
    static InputParameters validParams();

    ADStressBasedChemicalPotential(const InputParameters & parameters);
protected:

    ADReal computeQpResidual() override;
    
  /// Base name of the material system that this kernel applies to
  const std::string _base_name;

    /// The stress based chemical potential terms computed from the material for material swelling during intercalation
  const ADMaterialProperty<Real> & _swelling_chemical_potential;

    /// The stress based chemical potential terms computed from the material for the elastic energy term
  const ADMaterialProperty<Real> & _elastic_chemical_potential;


};