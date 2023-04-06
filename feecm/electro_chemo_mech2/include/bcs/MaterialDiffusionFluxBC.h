/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MaterialDiffusionFluxBC.h
 * Author: srinath
 *
 * Created on February 12, 2021, 5:22 PM
 */
#pragma once

#include "ADIntegratedBC.h"


class MaterialDiffusionFluxBC : public ADIntegratedBC
{
public:
    static InputParameters validParams();
    
    MaterialDiffusionFluxBC(const InputParameters & parameters);

protected:
    ADReal computeQpResidual() override;
    // Material prop diffusivity, TBD templated for values other than real

    const ADMaterialProperty<Real> * _diffusivity; 
};

