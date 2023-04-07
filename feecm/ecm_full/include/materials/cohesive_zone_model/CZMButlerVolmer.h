/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CZMButlerVolmer.h
 * Author: srinath
 *
 * Created on January 25, 2022, 4:04 PM
 */

#pragma once 
#include "CZMcomputeLocalFluxBase.h"

class CZMButlerVolmer : public CZMcomputeLocalFluxBase
{
    public: 
        static InputParameters validParams();
        CZMButlerVolmer(const InputParameters & parameters);
    protected:
        void computeInterfaceFluxAndDerivatives() override;
        const Real _d;
        const Function * _k_function;
        
        enum class Conductance
        {
          Conductance,
          Rct,
          Exchange_Current_density,
        } _conductance;        
        
        const bool _include_equil; 
        
        enum class SurfaceType
        {
          Primary,
          Secondary,
        } _surface_type;
        
        enum class Compute
        {
            Linear_Butler_Volmer,
            Butler_Volmer, 
            Lithium_Insertion,
        } _compute_type;
        
        const Real _cref; 
        const Real _faraday; 
        const Real _temperature; 
        const Real _gas_constant; 
        
        const Function * _reaction_rate_function; 
        MaterialProperty<Real> & _equilibrium_potential;
        const Real _RTF;
        const VariableValue * _coupled_value; 
        const Real _flux_sign;
};
