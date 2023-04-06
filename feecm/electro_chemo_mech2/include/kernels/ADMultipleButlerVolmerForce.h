/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   ADMultipleButlerVolmerForce.h
 * Author: srinath
 *
 * Created on March 24, 2022, 4:48 PM
 */

#pragma once

#include "ADKernel.h"
#include "ADButlerVolmerMaterial.h"

class ADMultipleButlerVolmerForce: public ADKernel
{
    public:
        static InputParameters validParams(); 
        ADMultipleButlerVolmerForce(const InputParameters & parameters);
    protected:
        virtual ADReal computeQpResidual() override;
//        const std::vector<std::string> _base_names;
        const Real _scale; // Scale the material property by this number
          /// number of materials
        const unsigned _num_materials;
        std::vector<const ADMaterialProperty<Real> *> _mat_properties;
//        std::vector<const ADButlerVolmerMaterial *> _materials;
};


