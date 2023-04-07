/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ADIstropicElasticSwelling.h
 * Author: srinath
 *
 * Created on January 16, 2021, 5:27 PM
 */
#pragma once

#include "ADIsoTropicHyperViscoSwellingBase.h"



class ADIsotropicElasticSwelling;


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

class ADIsotropicElasticSwelling : public ADIsoTropicHyperViscoSwellingBase
{
public:

	static InputParameters validParams();
  ADIsotropicElasticSwelling(const InputParameters & parameters);

protected:

    virtual void initQpStatefulProperties() override;
  
    virtual void propagateQpStatefulProperties() override;
    
    virtual void updateState(ADRankTwoTensor & strain_increment,
                           ADRankTwoTensor & inelastic_strain_increment,
                           const ADRankTwoTensor & rotation_increment,
                           ADRankTwoTensor & stress_new,
                           const RankTwoTensor & stress_old,
                           const ADRankFourTensor & elasticity_tensor,
                           const RankTwoTensor & elastic_strain_old) override;
//  virtual ADReal maximumPermissibleValue(const ADReal & effective_trial_stress) const override;
    
};
