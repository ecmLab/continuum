/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   SurfaceConcentration2ParameterAppox.h
 * Author: srinathcs
 *
 * Created on March 23, 2022, 9:15 AM
 */
#pragma once

#include "ADKernelValue.h"

class SurfaceConcentration2ParameterAppox: public ADKernelValue
{
    public:
        static InputParameters validParams();

        SurfaceConcentration2ParameterAppox(const InputParameters & parameters);
    protected:
        virtual ADReal precomputeQpResidual();
        const std::string _base_name;
        const ADVariableValue & _coupled_var;
    
        const ADMaterialProperty<Real> & _flux;
        const ADMaterialProperty<Real> & _diff;
        const ADMaterialProperty<Real> & _particle_size;
        const Real _scale;
};


