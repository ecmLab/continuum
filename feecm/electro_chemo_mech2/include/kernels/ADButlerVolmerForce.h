/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   ADButlerVolmerForce.h
 * Author: srinath
 *
 * Created on March 18, 2022, 9:07 AM
 */

#pragma once

#include "ADKernel.h"

class ADButlerVolmerForce: public ADKernel
{
    public:
        static InputParameters validParams(); 
        ADButlerVolmerForce(const InputParameters & parameters);
    protected:
        virtual ADReal computeQpResidual() override;
        const std::string _base_name;
        Real _scale; // Scale the material property by this number
        const ADMaterialProperty<Real> & _mat_property;
};

