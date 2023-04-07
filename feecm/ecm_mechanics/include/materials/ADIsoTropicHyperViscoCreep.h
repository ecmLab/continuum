
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADIsoTropicHyperViscoBase.h"

#define usingIsoTropicHyperViscoMembers                                                \
  usingIsoTropicHyperViscoBaseMembers


class ADIsoTropicHyperViscoCreep;


/**
 * This class uses the Discrete material in a radial return isotropic plasticity
 * model.  This class is one of the basic radial return constitutive models;
 * more complex constitutive models combine creep and plasticity.
 *
 * This class inherits from RadialReturnStressUpdate and must be used
 * in conjunction with ComputeReturnMappingStress.  This class calculates
 * an effective trial stress, an effective scalar plastic strain
 * increment, and the derivative of the scalar effective plastic strain increment;
 * these values are passed to the RadialReturnStressUpdate to compute
 * the radial return stress increment.  This isotropic plasticity class also
 * computes the plastic strain as a stateful material property.
 *
 * This class is based on the implicit integration algorithm in F. Dunne and N.
 * Petrinic's Introduction to Computational Plasticity (2004) Oxford University
 * Press, pg. 146 - 149.
 */

class ADIsoTropicHyperViscoCreep : public ADIsoTropicHyperViscoBase
{
public:
  static InputParameters validParams();

  ADIsoTropicHyperViscoCreep(const InputParameters & parameters);

protected:
    
    
//  virtual ADReal maximumPermissibleValue(const ADReal & effective_trial_stress) const override;
  
  virtual void initQpStatefulProperties() override;
  
  virtual void propagateQpStatefulProperties() override;
  
  virtual void computeStressInitialize(const ADReal & effective_trial_stress,
                                       const ADRankFourTensor & elasticity_tensor) override;
  
  virtual ADReal computeResidual(const ADReal & effective_trial_stress,
                                 const ADReal & scalar) override;
  
  virtual ADReal computeDerivative(const ADReal & effective_trial_stress,
                                   const ADReal & scalar) override;


  virtual void computeStressFinalize(const ADRankTwoTensor & plastic_strain_increment) override;
  virtual void iterationFinalize(ADReal scalar) override;

  virtual ADReal computeHardeningValue(const ADReal & scalar) ;
  virtual ADReal computeHardeningDerivative(const ADReal & scalar) ;
  
//   virtual Real computeReferenceResidual(const ADReal & effective_trial_stress,
//                                        const ADReal & scalar_effective_inelastic_strain) override;
   virtual Real computeTimeStepLimit() override;
 
  Real _mrate; // rate_exponent
  Real _Y0; // initial_resistance
  Real _H0; // Hardening modulus
  Real _Ysat; // Saturation value
  Real _ahard; // Hardening exponent

  Real _A;     // Pre-factor for activation energy
  Real _Q;     // Activation Energy
  Real _R; // Gas Constant
  Real _T;     // Temperature
  Real _n;     // Saturation hardening exponent
  
  Real _epsilon_rate; // reference_strain_rate
  
//  Real _alpha1, _alpha2, _alpha3, _omega;
  
  // Equivalent shear properties 
  Real _shear_rate;
  Real _shear_initial_resistance;
  Real _shear_saturation;
  Real _shear_initial_hardness;
  Real _prefactor;
  
  ADMaterialProperty<Real> & _saturation_strength;
  const MaterialProperty<Real> & _saturation_strength_old;

  
  template <typename T> 
  auto 
  sgn(T val) -> decltype(T(0)) 
  {
      return (T(0) < val) - (val < T(0));
  }
  
};
