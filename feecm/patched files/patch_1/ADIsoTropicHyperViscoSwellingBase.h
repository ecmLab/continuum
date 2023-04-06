
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "RadialReturnStressUpdate.h"

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
 * Directional Swelling model with
 * Large strain Logarithmic strain  based Istotropic Viscoplastic model
 * Based on Weber and Anand (1990) Computer Methods in App mech and engr
 * and Narayan and Anand (2020) Journal of Electrochemical Society"
 * Hardening a function of strain rate
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
                           const RankTwoTensor & elastic_strain_old, 
                           bool compute_full_tangent_operator,
                           RankFourTensor & tangent_operator) override;
  
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

   virtual void computeNormaltoBoundaryandPrincipalDirections();

   // Isotropic Swelling 
   bool _isotropic_swelling;

    // Rate form of swelling equations 
    bool _rate_form;
  
   // Concentration variable that induces swelling
   const ADVariableValue * _concentration;
   
   // Time derivative of concentration 
   const ADVariableValue * _cdot;
   
   // Reference concentration
   Real _cref;

   // Growth weights for swelling 
   RealVectorValue _alpha;
  
   // Molar volume of Li for swelling
   Real _omega;

   // Names of boundaries that the normals and principal directions are computed
   std::vector<BoundaryName> _boundary_name;
   
   // Flag for computing principal directions in deformed configuration 
   bool _use_deformed_directions;

  // a string to prepend to the plastic strain Material Property name
  const std::string _plastic_prepend;
  
  // plastic strain in this model
  ADMaterialProperty<RankTwoTensor> & _plastic_strain;

  // old value of plastic strain
  const MaterialProperty<RankTwoTensor> & _plastic_strain_old;
  
  // Hardening variable
  ADMaterialProperty<Real> & _hardening_variable;
  const MaterialProperty<Real> & _hardening_variable_old;
  
  // Value of current saturation strength
  ADMaterialProperty<Real> & _strength_variable;
  const MaterialProperty<Real> & _strength_variable_old;
  
  //  Deformation gradient inherited from FiniteStrain
  const ADMaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old;
  
  // The growth deformation gradient in FeFpFg decomposition
  ADMaterialProperty<RankTwoTensor> & _Fg;
  const MaterialProperty<RankTwoTensor> & _Fg_old;
  
  // Total volume change during swelling/growth
  ADMaterialProperty<Real> & _Jg;
  const MaterialProperty<Real> & _Jg_old;
  
  // Elastic portion of deformation gradient in FeFpFg decomposition
  ADMaterialProperty<RankTwoTensor> & _Fe;
  
  // Plastic portion of deformation gradient in FeFpFg decomposition
  ADMaterialProperty<RankTwoTensor> & _Fp;
  const MaterialProperty<RankTwoTensor> & _Fp_old;
  
  // Value of the current yield_strength
  ADMaterialProperty<Real> & _yield_strength;
  
  
  // The logarithmic elastic strain 
  ADMaterialProperty<RankTwoTensor> & _logarithmic_elastic_strain;
  const MaterialProperty<RankTwoTensor> & _logarithmic_elastic_strain_old;

  // The effective plastic strain rate
  ADMaterialProperty<Real> & _plastic_strain_rate;
  const MaterialProperty<Real> & _plastic_strain_rate_old;

  // The total effective plastic strain
  ADMaterialProperty<Real> & _effective_plastic_strain;
  const MaterialProperty<Real> & _effective_plastic_strain_old;

  // Principal direction for swelling/growth
  ADMaterialProperty<RealVectorValue> & _swell_normal;

  // Principal direction for swelling/growth
  ADMaterialProperty<RankTwoTensor> & _growth_tensor;

 
  // Value of the Elastic Mandel Stress 
  ADMaterialProperty<RankTwoTensor> & _mandel_stress;

  // Value of the Growth Stress (defined interms of Mandel stress)
  ADMaterialProperty<RankTwoTensor> & _growth_stress;

  // Value of the growth chemical potential
  ADMaterialProperty<Real> & _swelling_chemical_potential;

  // Value of the elastic chemical potential
  ADMaterialProperty<Real> & _elastic_chemical_potential;



  
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
  
//  const ADRankFourTensor _II;
//  
//  const ADRankFourTensor _II2;
  
};
