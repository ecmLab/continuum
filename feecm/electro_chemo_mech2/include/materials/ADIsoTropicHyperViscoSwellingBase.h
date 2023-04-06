    //* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADRadialReturnStressUpdate.h"


class ADIsoTropicHyperViscoSwellingBase;


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

class ADIsoTropicHyperViscoSwellingBase : public ADRadialReturnStressUpdate
{
public:
  static InputParameters validParams();

  ADIsoTropicHyperViscoSwellingBase(const InputParameters & parameters);

protected:
    
  virtual ADReal minimumPermissibleValue(const ADReal & effective_trial_stress) const override;
  
  virtual ADReal maximumPermissibleValue(const ADReal & effective_trial_stress) const override;

  
  
  virtual void updateState(ADRankTwoTensor & strain_increment,
                           ADRankTwoTensor & inelastic_strain_increment,
                           const ADRankTwoTensor & rotation_increment,
                           ADRankTwoTensor & stress_new,
                           const RankTwoTensor & stress_old,
                           const ADRankFourTensor & elasticity_tensor,
                           const RankTwoTensor & elastic_strain_old) override;
  
  virtual void initQpStatefulProperties() override ;
  virtual void propagateQpStatefulProperties() override;
  
  virtual ADReal initialGuess(const ADReal & /*effective_trial_stress*/) override 
  { return 1.0; }
  
  virtual void computeStressInitialize(const ADReal & effective_trial_stress,
                                       const ADRankFourTensor & elasticity_tensor) {};
  
  virtual ADReal computeResidual(const ADReal & effective_trial_stress,
                                 const ADReal & scalar) override { return 0; };
                                 
  virtual ADReal computeDerivative(const ADReal & effective_trial_stress,
                                   const ADReal & scalar) override {return 0; }; 
                                   
  virtual void iterationFinalize(ADReal scalar) {};
  virtual void computeStressFinalize(const ADRankTwoTensor & plastic_strain_increment) {};

  virtual Real computeTimeStepLimit() override;

   virtual void performKinematics(ADRankTwoTensor & R, ADRankTwoTensor  & E, 
                                  bool old_new);
     

   virtual void computeSwelling(ADReal & Jg, ADRankTwoTensor & Fg);
  
   
   const ADVariableValue * _concentration;
  
   Real _cref;

   Real _alpha1, _alpha2, _alpha3, _omega;
   
  /// a string to prepend to the plastic strain Material Property name
  const std::string _plastic_prepend;
  
   
  /// plastic strain in this model
  ADMaterialProperty<RankTwoTensor> & _plastic_strain;

  /// old value of plastic strain
  const MaterialProperty<RankTwoTensor> & _plastic_strain_old;
  
  
//
  ADMaterialProperty<Real> & _hardening_variable;
  const MaterialProperty<Real> & _hardening_variable_old;
//
  ADMaterialProperty<Real> & _strength_variable;
  const MaterialProperty<Real> & _strength_variable_old;
//  
  const ADMaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old;
  
  ADMaterialProperty<RankTwoTensor> & _Fg;
  ADMaterialProperty<Real> & _Jg;
  const MaterialProperty<Real> & _Jg_old;
//
  ADMaterialProperty<RankTwoTensor> & _Fp;
  const MaterialProperty<RankTwoTensor> & _Fp_old;
  
  ADMaterialProperty<Real> & _yield_strength;
  
  ADMaterialProperty<RankTwoTensor> & _logarithmic_elastic_strain;
  const MaterialProperty<RankTwoTensor> & _logarithmic_elastic_strain_old;

  ADMaterialProperty<Real> & _plastic_strain_rate;
  const MaterialProperty<Real> & _plastic_strain_rate_old;

  ADMaterialProperty<Real> & _effective_plastic_strain;
  const MaterialProperty<Real> & _effective_plastic_strain_old;
  
  ADMaterialProperty<RankTwoTensor> & _mandel_stress;

  
  /**
   * Rank two identity tensor
   */
  const ADRankTwoTensor _identity_two;

  /**
   * Rank four symmetric identity tensor
   */
  const ADRankFourTensor _identity_symmetric_four;

  /**
   * Rank four deviatoric projection tensor
   */
  const ADRankFourTensor _deviatoric_projection_four;
  
  const ADRankFourTensor _II;
  
  const ADRankFourTensor _II2;
  
};
