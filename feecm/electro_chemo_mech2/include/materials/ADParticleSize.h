/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   ADParticleSize.h
 * Author: srinathcs
 *
 * Created on March 23, 2022, 1:11 PM
 */
#pragma once
#include "Material.h"

class ADParticleSize: public ADMaterial
{
    public:
        static InputParameters validParams();

        ADParticleSize(const InputParameters & parameters);
//        Real surfacetoVolume();
    protected:
        virtual void computeQpProperties() override;
        const std::string _base_name;
        ADMaterialProperty<Real> & _particle_size;
        ADMaterialProperty<Real> & _surface_to_volume;
        const Real _size;
        const Function * _particle_size_function;
        enum class ParticleType
        {
            Spherical,
            Cylindrical,
        } _particle_type;
        const Real _volume_fraction;
        ADMaterialProperty<Real> & _vol_fraction;
};