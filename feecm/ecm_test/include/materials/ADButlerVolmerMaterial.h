/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   ADButlerVolmerMaterial.h
 * Author: srinath
 *
 * Created on March 18, 2022, 6:41 AM
 */

#pragma once 

#include "ADMaterial.h"

class ADButlerVolmerMaterial : public ADMaterial 
{
public:
    static InputParameters validParams();
    
    ADButlerVolmerMaterial(const InputParameters & parameters);
protected:
    virtual void computeQpProperties() override;
    const std::string _base_name;
    const Real _gas_constant;
    const Real _temperature; 
    const Real _faraday; 
    const ADVariableValue * _electrolyte_potential;
    const ADVariableValue * _electrode_potential;
    const bool _include_equil;
    const ADMaterialProperty<Real> * _equilibirium_potential;
    
    ADMaterialProperty<Real> & _bvcurrent; 
    ADMaterialProperty<Real> & _bvflux;
    ADMaterialProperty<Real> & _bvflux_force;
    ADMaterialProperty<Real> & _bvcurrent_force;
    ADMaterialProperty<Real> & _eta; 

    const ADMaterialProperty<Real> & _surface_to_volume; 
    const ADMaterialProperty<Real> & _volume_fraction;
       
    const ADMaterialProperty<Real> * _i0;
    const Real _a0; // Surface to volume ratio of cathode
private:

};


