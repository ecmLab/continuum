/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ADComputeEquilibriumPotential.h
 * Author: srinath
 *
 * Created on August 14, 2020, 12:43 PM
 */

#pragma once

#include "ADMaterial.h"
/**
 * Computes the equlibrium potential and chemical potential 
 * from the swelling_chemical_potential,
 * elastic_chemical_potential, concentration and reaction rate. 
 */


class ADComputeEquilibriumPotential:public ADMaterial {
public:
    static InputParameters validParams();
    
    ADComputeEquilibriumPotential(const InputParameters & parameters);
protected:
    virtual void computeQpProperties();
    
  /// Base name of the material system that this kernel applies to
  const std::string _base_name;
  ADReal _faraday; 
  ADReal _temp; 
  ADReal _gas_constant;
  bool _include_reaction_rate;
  const ADVariableValue * _concentration;
  bool _has_conc; 
  bool _include_conc; 
  Real _c0;
  bool _include_mechanical_effects;
  bool _exclude_elastic_energy;
  ADReal _reaction_rate;

  const Function * _reaction_rate_function; 
  bool _is_ocv;


    /// The stress based chemical potential terms computed from the material for material swelling during intercalation
  const ADMaterialProperty<Real> * const _swelling_chemical_potential;

    /// The stress based chemical potential terms computed from the material for the elastic energy term
  const ADMaterialProperty<Real> * const _elastic_chemical_potential;

  /// Reaction rate for nernst potential
  
   


  ADMaterialProperty<Real> & _equilibrium_potential;

  ADMaterialProperty<Real> & _chemical_potential;

  ADMaterialProperty<Real> & _state_of_charge;
  
    ///Diffusion potential  formulation
  enum class DiffPotential
  {
    DeHoff,
    Bower
  } _diff_potential;

private:

};


