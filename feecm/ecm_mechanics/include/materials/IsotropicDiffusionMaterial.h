/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   IsotropicDiffusionMaterial.h
 * Author: srinath
 *
 * Created on August 14, 2020, 12:43 PM
 */

#pragma once

#include "Material.h"
/**
 * IsotropicDiffusion Material
 * MaterialProperty that can be used as a Diffusion coefficient in diffusion 
 * calculations. It assumes that there are 3 orthonormal 
 * principal directions along which diffusion occurs. 
 */


class ADIsotropicDiffusionMaterial:public Material {
public:
    static InputParameters validParams();
    
    ADIsotropicDiffusionMaterial(const InputParameters & parameters);
protected:
   
   virtual void computeQpProperties();

  ADReal _diffusion_coef;
   /// Diffusivity tensor material property
  ADMaterialProperty<Real> & _diffusivity;

private:

};


