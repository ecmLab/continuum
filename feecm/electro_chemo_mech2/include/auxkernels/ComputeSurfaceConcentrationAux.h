/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   ComputeSurfaceConcentrationAux.h
 * Author: srinath
 *
 * Created on March 21, 2022, 9:09 AM
 */

#pragma once 
#include "AuxKernel.h"

class ComputeSurfaceConcentrationAux;

class ComputeSurfaceConcentrationAux : public AuxKernel 
{
    public:
        static InputParameters validParams();
    
        ComputeSurfaceConcentrationAux(const InputParameters & parameters);

    protected:
        virtual Real computeValue();
        const ADVariableValue & _concentration; 
        const ADMaterialProperty<Real> & _flux; 
        const Real _particle_size; 
        const ADMaterialProperty<Real> & _diffusion_coefficient;
        const Real _scale; 
};

