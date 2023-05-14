/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   ADCoupledButlerVolmerForce.h
 * Author: srinathcs
 *
 * Created on March 17, 2022, 11:33 AM
 */
#pragma once
#include "ADKernel.h"

class ADCoupledButlerVolmerForce : public ADKernel

{
    public:
        static InputParameters validParams();
        
        ADCoupledButlerVolmerForce(const InputParameters & parameters);
    protected:
        virtual ADReal computeQpResidual() override;
    private:
        Real _i0;    // exchange (equilibrium) reaction rate
        Real _faraday;    // Faraday constant (96485.33289 C/mol)
        Real _gas_constant;    // gas constant (8.3144598 J/mol/K)
        Real _temperature;    // temperature
//        Real _coefficient;    // concentration/coefficent

        // Coupled variable number of electrolyte potential
        unsigned int _num_electrolyte_potential_var;
        const ADVariableValue & _electrolyte_potential;
        // Coupled variable number of electrolyte potential
        unsigned int _num_electrode_potential_var;
        const ADVariableValue & _electrode_potential;
        const ADMaterialProperty<Real> & _equilibrium_potential;
        
};