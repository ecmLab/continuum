/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   GaoBower3DCTraction.h
 * Author: srinath
 *
 * Created on January 27, 2022, 12:09 PM
 */

#pragma once

#include "SalehaniIrani3DCTraction.h"

#include "GaoBower3DCTraction.h"

class GaoBower3DCTraction : public SalehaniIrani3DCTraction
{
    public:
        static InputParameters validParams();
        GaoBower3DCTraction(const InputParameters & parameters);
    protected:
         void computeInterfaceTractionAndDerivatives() override;

        /// method computing the total traction
        RealVectorValue computeTraction();

        /// method computing the total traction derivatives w.r.t. the interface displacement jump
        RankTwoTensor computeTractionDerivatives();
        
        /// Normal and tangential viscosity
        
        const RealVectorValue _viscosity;
        
        /// Xu-Needleman parameters
        const Real _r; 
        const Real _q; 
        
        /// The displacment jump  incremenet in local coordinates
        MaterialProperty<RealVectorValue> & _interface_displacement_jump_inc;

        /// The old interface displacment jump
        const MaterialProperty<RealVectorValue> & _interface_displacement_jump_old;

};