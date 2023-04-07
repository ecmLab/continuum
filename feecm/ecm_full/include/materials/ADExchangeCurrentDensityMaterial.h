/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   ADExchangeCurrentDensityMaterial.h
 * Author: srinath
 *
 * Created on March 21, 2022, 3:51 PM
 */

#pragma once

#include "Material.h"

class ADExchangeCurrentDensityMaterial : public ADMaterial
{
    public:
        static InputParameters validParams();
    
        ADExchangeCurrentDensityMaterial(const InputParameters & parameters);
    protected:
        virtual void computeQpProperties() override;
        const std::string _base_name;
        ADMaterialProperty<Real> & _i0;
        enum class ExchangeCurrentDensityType
        {
            Constant,
            Lithium_Insertion
        } _exchange_type;
        const ADVariableValue & _concentration; 
        const Real _i_ref;
        const Real _cmax;
};