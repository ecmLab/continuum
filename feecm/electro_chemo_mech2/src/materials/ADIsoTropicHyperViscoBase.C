//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADIsoTropicHyperViscoBase.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

#include "ElasticityTensorTools.h"

//registerADMooseObject("electro_chemo_mechApp", ADIsoTropicHyperViscoBase);

InputParameters
ADIsoTropicHyperViscoBase::validParams()
{
	InputParameters params = ADRadialReturnStressUpdate::validParams();
	params.addClassDescription("Large strain Logarithmic strain based Istotropic Viscoplastic model "
	                               "Based on Weber and Anand (1990) Computer Methods in App mech and engr"
	                               "Hardening a function of strain rate");
	        params.addDeprecatedParam<std::string>(
	        "plastic_prepend",
	        "",
	        "String that is prepended to the plastic_strain Material Property",
	        "This has been replaced by the 'base_name' parameter");
	params.set<std::string>("effective_inelastic_strain_name") = "effective_plastic_strain";
	return params;
}

ADIsoTropicHyperViscoBase::ADIsoTropicHyperViscoBase(
    const InputParameters & parameters)
  : ADRadialReturnStressUpdate(parameters),
    _plastic_prepend(getParam<std::string>("plastic_prepend")),
    _plastic_strain(
        declareADProperty<RankTwoTensor>(_base_name + _plastic_prepend + "plastic_strain")),
    _plastic_strain_old(
        getMaterialPropertyOld<RankTwoTensor>(_base_name + _plastic_prepend + "plastic_strain")),
    _hardening_variable(declareADProperty<Real>(_base_name + "hardening_variable")),
    _hardening_variable_old(getMaterialPropertyOld<Real>(_base_name + "hardening_variable")),
    _strength_variable(declareADProperty<Real>(_base_name + "strength_variable")),
    _strength_variable_old(getMaterialPropertyOld<Real>(_base_name + "strength_variable")),
    _deformation_gradient(getADMaterialProperty<RankTwoTensor>(_base_name + "deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "deformation_gradient")),
    _Fp(declareADProperty<RankTwoTensor>(_base_name + _plastic_prepend + "plastic_distortion")),
    _Fp_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + _plastic_prepend + "plastic_distortion")),
    _yield_strength(declareADProperty<Real>(_base_name + "yield_strength")), 
    _logarithmic_elastic_strain(declareADProperty<RankTwoTensor>(_base_name + _plastic_prepend + "logarithmic_elastic_strain")),
    _logarithmic_elastic_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + _plastic_prepend + "logarithmic_elastic_strain")),
    _plastic_strain_rate(declareADProperty<Real>(_base_name + "plastic_strain_rate")),
    _plastic_strain_rate_old(getMaterialPropertyOld<Real>(_base_name + "plastic_strain_rate")),
    _effective_plastic_strain(declareADProperty<Real>(_base_name + "effective_plastic_strain")),
    _effective_plastic_strain_old(getMaterialPropertyOld<Real>(_base_name + "effective_plastic_strain")),
    _mandel_stress(declareADProperty<RankTwoTensor>(_base_name + "mandel_stress")),        
    _identity_two(ADRankTwoTensor::initIdentity),
    _identity_symmetric_four(ADRankFourTensor::initIdentitySymmetricFour),
    _deviatoric_projection_four(_identity_symmetric_four -
                                _identity_two.outerProduct(_identity_two) / 3.0),
    _II(_identity_two.mixedProductIkJl(_identity_two) 
                          + _identity_two.mixedProductIlJk(_identity_two)),
    _II2(_identity_two.outerProduct(_identity_two))
        
{
//    _mandel_stress[_qp].zero();
//    _Fp[_qp].zero();
//    _Fp[_qp].addIa(1.0);
}


ADReal
ADIsoTropicHyperViscoBase::minimumPermissibleValue(
    const ADReal & /*effective_trial_stress*/) const
{
  return 0.0;
}

ADReal
ADIsoTropicHyperViscoBase::maximumPermissibleValue(
    const ADReal & effective_trial_stress) const
{ 
  return 10000.0;
}

void
ADIsoTropicHyperViscoBase::initQpStatefulProperties()
{
  _plastic_strain[_qp].zero();
  _Fp[_qp].zero();
  _Fp[_qp].addIa(1.0);
  
  _logarithmic_elastic_strain[_qp].zero();
  _plastic_strain_rate[_qp] = 0.0;
  _effective_plastic_strain[_qp] = 0.0;
}

void
ADIsoTropicHyperViscoBase::propagateQpStatefulProperties()
{
  _plastic_strain[_qp] = _plastic_strain_old[_qp];
  _Fp[_qp] = _Fp_old[_qp];  
  _logarithmic_elastic_strain[_qp] = _logarithmic_elastic_strain_old[_qp];
  _plastic_strain_rate[_qp] = _plastic_strain_rate_old[_qp];
  _effective_plastic_strain[_qp]  = _effective_plastic_strain_old[_qp];
  propagateQpStatefulPropertiesRadialReturn();
}


void
ADIsoTropicHyperViscoBase::updateState(
    ADRankTwoTensor & strain_increment,
    ADRankTwoTensor & inelastic_strain_increment,
    const ADRankTwoTensor & /*rotation_increment*/,
    ADRankTwoTensor & stress_new,
    const RankTwoTensor & /*stress_old*/,
    const ADRankFourTensor & elasticity_tensor,
    const RankTwoTensor & elastic_strain_old)
{
    
   // compute the deviatoric trial stress and trial strain from the current intermediate
  // configuration
   
    const ADReal mu = ElasticityTensorTools::getIsotropicShearModulus(elasticity_tensor);
    const ADReal bulk  = ElasticityTensorTools::getIsotropicBulkModulus(elasticity_tensor);
          
    ADRankTwoTensor Ee_tr, R_tr, Fp;
//    Fp = _Fp_old[_qp];
    performKinematics(R_tr, Ee_tr, false);

    // Mandel Stress trial
    ADRankTwoTensor Me_tr = (mu * _II + (bulk - 2.0/3.0 * mu ) * _II2) * Ee_tr;
    ADRankTwoTensor deviatoric_trial_stress = Me_tr.deviatoric();   
    
//  // compute the effective trial stress
  ADReal dev_trial_stress_squared =
      deviatoric_trial_stress.doubleContraction(deviatoric_trial_stress);
  
  ADReal norm_dev_trial_stress = MooseUtils::absoluteFuzzyEqual(dev_trial_stress_squared, 0.0)
                                      ? 0.0
                                      : std::sqrt(dev_trial_stress_squared);
  ADReal effective_trial_stress = MooseUtils::absoluteFuzzyEqual(dev_trial_stress_squared, 0.0)
                                      ? 0.0
                                      : std::sqrt(1.0 / 2.0 * dev_trial_stress_squared);


    ADRankTwoTensor flow_direction(ADRankTwoTensor::initIdentity);

  // Set the value of 3 * shear modulus for use as a reference residual value
  _three_shear_modulus = 3.0 * ElasticityTensorTools::getIsotropicShearModulus(elasticity_tensor);

  computeStressInitialize(effective_trial_stress, elasticity_tensor);

  // Use Newton iteration to determine the scalar effective inelastic strain increment
  ADReal scalar_effective_inelastic_strain = 0.0;
  if (!MooseUtils::absoluteFuzzyEqual(effective_trial_stress, 0.0))
  {
    flow_direction = deviatoric_trial_stress / (std::sqrt(2.0) * effective_trial_stress);
    returnMappingSolve(effective_trial_stress, scalar_effective_inelastic_strain, _console);
    if (scalar_effective_inelastic_strain != 0.0)
      inelastic_strain_increment = std::sqrt(1.0/2.0) * scalar_effective_inelastic_strain * 
              flow_direction;
    else
      inelastic_strain_increment.zero();
  }
  else
    inelastic_strain_increment.zero();
  ADRankTwoTensor dpt = inelastic_strain_increment * _dt;


    strain_increment -= inelastic_strain_increment;
  // Use the old elastic strain here because we require tensors used by this class
  // to be isotropic and this method natively allows for changing in time
  // elasticity tensors
  if (scalar_effective_inelastic_strain > 0.0) 
  {
    _plastic_strain_rate[_qp] = scalar_effective_inelastic_strain;
    _effective_inelastic_strain[_qp] =
      _effective_inelastic_strain_old[_qp] + _dt * 
          scalar_effective_inelastic_strain / std::sqrt(3.0);
    _effective_plastic_strain[_qp]  = _effective_plastic_strain_old[_qp] +   
            _dt * scalar_effective_inelastic_strain / std::sqrt(3.0);
  //  // Compute new Cauchy stress
    std::vector<ADReal> e_value(3);
    ADRankTwoTensor e_vector, N1, N2, N3;

    dpt.symmetricEigenvaluesEigenvectors(e_value, e_vector);

    N1.vectorOuterProduct(e_vector.column(0), e_vector.column(0));
    N2.vectorOuterProduct(e_vector.column(1), e_vector.column(1));
    N3.vectorOuterProduct(e_vector.column(2), e_vector.column(2));
    ADRankTwoTensor expDp = N1 * std::exp(e_value[0]) + 
                          N2 * std::exp(e_value[1]) + 
                          N3 * std::exp(e_value[2]);  
    _Fp[_qp] = expDp * _Fp_old[_qp];

    mooseAssert(_Fp[_qp].det() > 0.0, "Det of Fp less than 0");
    
  } else
  {
      _plastic_strain_rate[_qp]  = _plastic_strain_rate_old[_qp];
      _effective_inelastic_strain[_qp] = _effective_inelastic_strain_old[_qp];
      _effective_plastic_strain[_qp]  = _effective_plastic_strain_old[_qp];
      _Fp[_qp] = _Fp_old[_qp];
  }
    R_tr.zero();
    Ee_tr.zero();
    Fp = _Fp[_qp];
  performKinematics(R_tr, Ee_tr, true);
  ADRankTwoTensor Me_tau = Me_tr - (2.0/3.0)*_three_shear_modulus 
                         * _dt * inelastic_strain_increment;
  
  _logarithmic_elastic_strain[_qp] = Ee_tr;
  _mandel_stress[_qp] = Me_tau;
  stress_new = R_tr * (Me_tau * R_tr.transpose()) / _deformation_gradient[_qp].det();
  
  
//  stress_new = elasticity_tensor * (elastic_strain_old + strain_increment);

  computeStressFinalize(inelastic_strain_increment);
}


Real
ADIsoTropicHyperViscoBase:: computeTimeStepLimit()
{
    Real scalar_inelastic_strain_incr;
    
    scalar_inelastic_strain_incr = std::abs (MetaPhysicL::raw_value(_effective_inelastic_strain[_qp] -
            _effective_inelastic_strain_old[_qp]))/ _max_inelastic_increment;
    if (MooseUtils::absoluteFuzzyEqual(scalar_inelastic_strain_incr, 0.0))
        return std::numeric_limits<Real>::max();

    if (scalar_inelastic_strain_incr <= 0.5)
        return _dt * 1.5;
    else if (scalar_inelastic_strain_incr > 0.5 && scalar_inelastic_strain_incr <= 0.8)
        return _dt * 1.25;
    else if (scalar_inelastic_strain_incr > 0.8 && scalar_inelastic_strain_incr <= 1.25)
        return _dt * 0.75;
    else
        return _dt * 0.5;
}

void
ADIsoTropicHyperViscoBase::performKinematics(ADRankTwoTensor & R,
        ADRankTwoTensor & E, bool old_new)
{
    ADRankTwoTensor Fe;
    ADRankTwoTensor Fi;
    if (old_new)
    {
        Fi = _Fp[_qp];
    }
    else
    {
        Fi = _Fp_old[_qp];
    }
    
    
//    if ( Fi.det() <= 0)
//    {
//        Fi.zero();
//        Fi.addIa(1.0);
//    }
    Fe = _deformation_gradient[_qp]*Fi.inverse();
    
    ADRankTwoTensor Ce = Fe.transpose() * Fe; 
    std::vector<ADReal> e_value(3);
    ADRankTwoTensor e_vector, N1, N2, N3;
    
    Ce.symmetricEigenvaluesEigenvectors(e_value, e_vector);
    
    const auto lambda1 = std::sqrt(e_value[0]);
    const auto lambda2 = std::sqrt(e_value[1]);
    const auto lambda3 = std::sqrt(e_value[2]);
    
    N1.vectorOuterProduct(e_vector.column(0), e_vector.column(0));
    N2.vectorOuterProduct(e_vector.column(1), e_vector.column(1));
    N3.vectorOuterProduct(e_vector.column(2), e_vector.column(2));
    
    ADRankTwoTensor U = N1 * lambda1 + N2 * lambda2 + N3 * lambda3;
    ADRankTwoTensor invUe(U.inverse());
    E = N1 * std::log(lambda1) + N2 * std::log(lambda2) + N3 * std::log(lambda3);
    R = Fe * invUe;

}

