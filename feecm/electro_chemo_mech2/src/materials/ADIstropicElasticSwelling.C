/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ADIstropicElasticSwelling.C
 * Author: srinath
 * 
 * Created on January 16, 2021, 5:27 PM
 */

#include "ADIstropicElasticSwelling.h"

#include "ElasticityTensorTools.h"

registerADMooseObject("electro_chemo_mechApp", ADIsotropicElasticSwelling);


InputParameters
ADIsotropicElasticSwelling::validParams()
{
	InputParameters params = ADIsoTropicHyperViscoSwellingBase::validParams();
    params.addClassDescription("Large strain Logarithmic strain based Istotropic Viscoplastic model "
                               "Based on Weber and Anand (1990) Computer Methods in App mech and engr"
                               "Hardening a function of strain rate");
    params.set<bool>("elastic_material") = true;
    params.set<std::string>("effective_inelastic_strain_name") = "effective_plastic_strain";
    return params;
}


ADIsotropicElasticSwelling::ADIsotropicElasticSwelling(const InputParameters & parameters)
  : ADIsoTropicHyperViscoSwellingBase(parameters)

{
}

void ADIsotropicElasticSwelling::initQpStatefulProperties()
{
      ADIsoTropicHyperViscoSwellingBase::initQpStatefulProperties();
}

void
ADIsotropicElasticSwelling::propagateQpStatefulProperties()
{
  ADIsoTropicHyperViscoSwellingBase::propagateQpStatefulProperties();
}


void
ADIsotropicElasticSwelling::updateState(
    ADRankTwoTensor & strain_increment,
    ADRankTwoTensor & inelastic_strain_increment,
    const ADRankTwoTensor & /*rotation_increment*/,
    ADRankTwoTensor & stress_new,
    const RankTwoTensor & /*stress_old*/,
    const ADRankFourTensor & elasticity_tensor,
    const RankTwoTensor & elastic_strain_old, 
    bool /*compute_full_tangent_operator = false*/,
    RankFourTensor & /*tangent_operator = _identityTensor*/)
{
    const ADReal mu = ElasticityTensorTools::getIsotropicShearModulus(elasticity_tensor);
    const ADReal bulk  = ElasticityTensorTools::getIsotropicBulkModulus(elasticity_tensor);
// compute the swelling deformation gradient
    ADReal Jg; 
    ADRankTwoTensor Fg; 
    Jg = 0.0;
    Fg.zero();
    ADIsoTropicHyperViscoSwellingBase::computeSwelling(Jg, Fg);
    _Jg[_qp] = Jg;
    _Fg[_qp] = Fg;
// compute the elastic deformation gradient
    ADRankTwoTensor Ee_tr, R_tr;
    ADIsoTropicHyperViscoSwellingBase::performKinematics(R_tr, Ee_tr, false);
//   // Trial Mandel Stress
//    ADReal trE0 = std::log(_deformation_gradient[_qp].det()/_Jg[_qp]);
//    // This has to be done because the trace of the logarithmic elasticity tensor 
//    // is the new one. 
//    ADRankTwoTensor Me_tr = 2.0 * mu * Ee_tr.deviatoric() + bulk * trE0 * _identity_two;
  ADRankTwoTensor Me_tau = (mu * _II + (bulk - 2.0/3.0 * mu ) * _II2) * Ee_tr;  
  
  _logarithmic_elastic_strain[_qp] = Ee_tr;
  _mandel_stress[_qp] = Me_tau;
  
  stress_new = R_tr * (Me_tau * R_tr.transpose()) 
          / _deformation_gradient[_qp].det() * _Jg[_qp];
   

//  stress_new = elasticity_tensor * (elastic_strain_old + strain_increment);
  _growth_stress[_qp] = _mandel_stress[_qp];
  /// We cannot apply any mechanical stresses on a material that does not exist
  if ((*_concentration)[_qp] > 0){
    _swelling_chemical_potential[_qp] = -_omega * _growth_stress[_qp].doubleContraction(_growth_tensor[_qp]);
    _elastic_chemical_potential[_qp] = _omega * 0.5 * _logarithmic_elastic_strain[_qp].doubleContraction(
            elasticity_tensor * _logarithmic_elastic_strain[_qp]); 
  } else
  {
      _swelling_chemical_potential[_qp]= 0.0;
      _elastic_chemical_potential[_qp] = 0.0;
              
  }

}