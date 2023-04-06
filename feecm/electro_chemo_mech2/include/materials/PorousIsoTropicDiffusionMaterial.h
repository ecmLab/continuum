/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   PorousIsoTropicDiffusionMaterial.h
 * Author: srinathcs
 *
 * Created on March 23, 2022, 4:25 PM
 */
#pragma once
#include "IsotropicDiffusionMaterial.h"

class PorousIsotropicDiffusionMaterial : public ADIsotropicDiffusionMaterial
{   
    public:
        static InputParameters validParams(); 
        PorousIsotropicDiffusionMaterial(const InputParameters & parameters);
    protected:
        virtual Real computeQpCorrection() override;
        const Real _volume_fraction;
        const Real _bruggeman_factor;
        enum class DiffusionModel
        {
           Bruggeman,
           Tortuosity, 
           User_Defined,
        } _diffusion_model;
        const Real _tortuosity;
        const Real _correction_factor;
    
};
