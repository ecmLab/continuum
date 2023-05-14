/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   DiffusionAlongPrincipalDirectionsMaterial.h
 * Author: srinath
 *
 * Created on August 14, 2020, 12:43 PM
 */

#pragma once

#include "Material.h"
/**
 * DiffusionAlongPrincipalDirections provides a simple RealTensorValue type 
 * MaterialProperty that can be used as a Diffusion coefficient in diffusion 
 * calculations. It assumes that there are 3 orthonormal 
 * principal directions along which diffusion occurs. 
 */


class ADDiffusionAlongPrincipalDirectionsMaterial:public Material {
public:
    static InputParameters validParams();
    
    ADDiffusionAlongPrincipalDirectionsMaterial(const InputParameters & parameters);
protected:
    virtual void computeQpProperties();
   const std::string _base_name; 
   // Mobility/diffusivity values along principal directions
   RealVectorValue _diffusivity_values;
   
  /// Material Property of the normal direction to swelling
  const ADMaterialProperty<RealVectorValue> & _swell_normal;

  
   /// Diffusivity tensor material property
  ADMaterialProperty<RealTensorValue> & _diffusivity;

private:

};


