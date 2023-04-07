//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADIsoTropicHyperViscoSwellingBase.h"
#include "libmesh/quadrature.h"
#include "ElasticityTensorTools.h"
#include "MooseMesh.h"

//registerADMooseObject("electro_chemo_mechApp", ADIsoTropicHyperViscoSwellingBase);

InputParameters
ADIsoTropicHyperViscoSwellingBase::validParams()
{
	InputParameters params = ADRadialReturnStressUpdate::validParams();

    params.addClassDescription("Directional Swelling model with "
                               "Large strain Logarithmic strain  based Istotropic Viscoplastic model "
                               "Based on Weber and Anand (1990) Computer Methods in App mech and engr"
                               "and Narayan and Anand (2020) Journal of Electrochemical Society"
                               "Hardening a function of strain rate");
        params.addParam<bool>("isotropic_swelling",false, "Isotropic Swelling");
        params.addParam<bool>("rate_form", false, "Rate form of Swelling equations");
        params.addCoupledVar("concentration", "Coupled concentration");
        params.addParam<Real>("cref",0, "Reference concentration");
        params.addRequiredParam<RealVectorValue>("alpha", "alpha");        
        params.addRequiredParam<Real>("omega", "Molar Volume of Li");
        params.addParam<std::vector<BoundaryName>>(
          "intBnd", "The boundary names of the interface. If no boundary is given"
                ",the default normal direction is assumed (0 1 0) global y-axis");
        params.addParam<bool>("use_deformed_directions",false, "Compute principal" ""
                                    "directions in deformed configuration");
        params.addDeprecatedParam<std::string>(
        "plastic_prepend",
        "",
        "String that is prepended to the plastic_strain Material Property",
        "This has been replaced by the 'base_name' parameter");
    params.set<std::string>("effective_inelastic_strain_name") = "effective_plastic_strain";
    return params;
}

ADIsoTropicHyperViscoSwellingBase::ADIsoTropicHyperViscoSwellingBase(const InputParameters & parameters)
  : ADRadialReturnStressUpdate(parameters),
    _isotropic_swelling(getParam<bool>("isotropic_swelling")),
    _rate_form(getParam<bool>("rate_form")),
    _concentration(isParamValid("concentration") ? &adCoupledValue("concentration") : nullptr),
    _cdot(isParamValid("concentration") ? &adCoupledDot("concentration") : nullptr),
    _cref(getParam<Real>("cref")),
    _alpha(getParam<RealVectorValue>("alpha")),
    _omega(getParam<Real>("omega")),
    _boundary_name(getParam<std::vector<BoundaryName>>("intBnd")),
    _use_deformed_directions(getParam<bool>("use_deformed_directions")),
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
    _Fg(declareADProperty<RankTwoTensor>(_base_name + _plastic_prepend + "growth_deformation")),
    _Fg_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + _plastic_prepend + "growth_deformation")),
    _Jg(declareADProperty<Real>(_base_name + _plastic_prepend + "swelling_vol_change")),    
    _Jg_old(getMaterialPropertyOld<Real>(_base_name + _plastic_prepend + "swelling_vol_change")),        
    _Fe(declareADProperty<RankTwoTensor>(_base_name + _plastic_prepend + "elastic_distortion")),
    _Fp(declareADProperty<RankTwoTensor>(_base_name + _plastic_prepend + "plastic_distortion")),
    _Fp_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + _plastic_prepend + "plastic_distortion")),
    _yield_strength(declareADProperty<Real>(_base_name + "yield_strength")), 
    _logarithmic_elastic_strain(declareADProperty<RankTwoTensor>(_base_name + _plastic_prepend + "logarithmic_elastic_strain")),
    _logarithmic_elastic_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + _plastic_prepend + "logarithmic_elastic_strain")),
    _plastic_strain_rate(declareADProperty<Real>(_base_name + "plastic_strain_rate")),
    _plastic_strain_rate_old(getMaterialPropertyOld<Real>(_base_name + "plastic_strain_rate")),
    _effective_plastic_strain(declareADProperty<Real>(_base_name + "effective_plastic_strain")),
    _effective_plastic_strain_old(getMaterialPropertyOld<Real>(_base_name + "effective_plastic_strain")),
    _swell_normal(declareADProperty<RealVectorValue>(_base_name + "swell_normal")),        
    _growth_tensor(declareADProperty<RankTwoTensor>(_base_name + "growth_tensor")),
    _mandel_stress(declareADProperty<RankTwoTensor>(_base_name + "mandel_stress")),
    _growth_stress(declareADProperty<RankTwoTensor>(_base_name + "growth_stress")),        
    _swelling_chemical_potential(declareADProperty<Real>(_base_name + "swelling_chemical_potential")),        
    _elastic_chemical_potential(declareADProperty<Real>(_base_name + "elastic_chemical_potential")),        
    _identity_two(ADRankTwoTensor::initIdentity),
    _identity_symmetric_four(ADRankFourTensor::initIdentitySymmetricFour),
    _deviatoric_projection_four(_identity_symmetric_four -
                                _identity_two.outerProduct(_identity_two) / 3.0)
//    _II(_identity_two.mixedProductIkJl(_identity_two) 
//                          + _identity_two.mixedProductIlJk(_identity_two)),
//    _II2(_identity_two.outerProduct(_identity_two))
        
{
    if (_isotropic_swelling)
    {
            _alpha(0) = 1.0/3.0;
            _alpha(1) = 1.0/3.0;
            _alpha(2) = 1.0/3.0;
    }
    ADReal sum_alpha = _alpha(0) + _alpha(1) + _alpha(2);
//    if (!_isotropic_swelling)
//    {
        if(!MooseUtils::absoluteFuzzyEqual(sum_alpha, 1.0, 1.0e-3))
            mooseError("Sum of swelling exponents has to be equal to 1.0");
//    }
   
}

ADReal
ADIsoTropicHyperViscoSwellingBase::minimumPermissibleValue(const ADReal & /*effective_trial_stress*/) const
{
  return 0.0;
}

ADReal
ADIsoTropicHyperViscoSwellingBase::maximumPermissibleValue(const ADReal & /*effective_trial_stress*/) const
{
  return 10000.0;
}


void
ADIsoTropicHyperViscoSwellingBase::initQpStatefulProperties()
{
  _plastic_strain[_qp].zero();
  _Fe[_qp].zero();
  _Fe[_qp].addIa(1.0);
  _Fp[_qp].zero();
  _Fp[_qp].addIa(1.0);
  _Fg[_qp].zero();
  _Fg[_qp].addIa(1.0);
  _Jg[_qp] = 1.0;
  
  _logarithmic_elastic_strain[_qp].zero();
  _plastic_strain_rate[_qp] = 0.0;
  _effective_plastic_strain[_qp] = 0.0;

// Compute normals of _qp points for swelling
  if (!_isotropic_swelling)
    computeNormaltoBoundaryandPrincipalDirections();

}

void
ADIsoTropicHyperViscoSwellingBase::propagateQpStatefulProperties()
{
  _plastic_strain[_qp] = _plastic_strain_old[_qp];
  _Fp[_qp] = _Fp_old[_qp];  
  _Jg[_qp] = _Jg_old[_qp];
  _Fg[_qp] = _Fg_old[_qp];
  _logarithmic_elastic_strain[_qp] = _logarithmic_elastic_strain_old[_qp];
  _plastic_strain_rate[_qp] = _plastic_strain_rate_old[_qp];
  _effective_plastic_strain[_qp]  = _effective_plastic_strain_old[_qp];
  propagateQpStatefulPropertiesRadialReturn();
}


void
ADIsoTropicHyperViscoSwellingBase::updateState(
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
// compute the deviatoric trial stress and trial strain from the current intermediate
// configuration
   
    const ADReal mu = ElasticityTensorTools::getIsotropicShearModulus(elasticity_tensor);
    const ADReal bulk  = ElasticityTensorTools::getIsotropicBulkModulus(elasticity_tensor);
// compute the swelling deformation gradient
    ADReal Jg; 
    ADRankTwoTensor Fg; 
    Jg = 0.0;
    Fg.zero();
    computeSwelling(Jg, Fg);
    _Jg[_qp] = Jg;
    _Fg[_qp] = Fg;
// compute the elastic deformation gradient
    ADRankTwoTensor Ee_tr, R_tr;
    performKinematics(R_tr, Ee_tr, false);

   // Trial Mandel Stress
    ADReal trE0 = std::log(_deformation_gradient[_qp].det()/_Jg[_qp]);
    // This has to be done because the trace of the logarithmic elasticity tensor 
    // is the new one. 
    ADRankTwoTensor Me_tr = 2.0 * mu * Ee_tr.deviatoric() + bulk * trE0 * _identity_two;

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

  Real eff_tr_stress = MetaPhysicL::raw_value(effective_trial_stress);
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
  performKinematics(R_tr, Ee_tr, true );
  
//  ADRankTwoTensor Me_tau = (mu * _II + (bulk - 2.0/3.0 * mu ) * _II2) * Ee_tr;
  ADRankTwoTensor Me_tau = Me_tr - (2.0/3.0)*_three_shear_modulus 
                         * _dt * inelastic_strain_increment;
  
  
  _logarithmic_elastic_strain[_qp] = Ee_tr;
  _mandel_stress[_qp] = Me_tau;
  
  stress_new = R_tr * (Me_tau * R_tr.transpose()) 
          / _deformation_gradient[_qp].det() * _Jg[_qp];
   

//  stress_new = elasticity_tensor * (elastic_strain_old + strain_increment);
   _growth_stress[_qp] = _Fp[_qp].transpose() * _mandel_stress[_qp] * (_Fp[_qp].inverse()).transpose();
  if ((*_concentration)[_qp] > 0)
  {
      _swelling_chemical_potential[_qp] = -_omega * _growth_stress[_qp].doubleContraction(_growth_tensor[_qp]);
    _elastic_chemical_potential[_qp] = _omega * 0.5 * _logarithmic_elastic_strain[_qp].doubleContraction(
               elasticity_tensor * _logarithmic_elastic_strain[_qp]); 
  } else
  {
      _swelling_chemical_potential[_qp] = 0.0;
      _elastic_chemical_potential[_qp] = 0.0;
  }
   
  computeStressFinalize(inelastic_strain_increment);
}



Real
ADIsoTropicHyperViscoSwellingBase:: computeTimeStepLimit()
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
ADIsoTropicHyperViscoSwellingBase::performKinematics(ADRankTwoTensor & R,
        ADRankTwoTensor & E, bool old_new)
{
//    ADRankTwoTensor Fe;
    ADRankTwoTensor Fi;
    if (old_new)
    {
        Fi = _Fp[_qp] * _Fg[_qp];
    }
    else
    {
        Fi = _Fp_old[_qp] * _Fg[_qp];
    }
    
    _Fe[_qp] = _deformation_gradient[_qp]*Fi.inverse();
    
    ADRankTwoTensor Ce = _Fe[_qp].transpose() * _Fe[_qp]; 
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
    R = _Fe[_qp] * invUe;

}


void
ADIsoTropicHyperViscoSwellingBase::computeSwelling(
                                    ADReal & Jg, ADRankTwoTensor & Fg)
{
    // This assumes that the 3 principal directions are along the cartesian plane
      // Any 3 orthogonal vectors can be defined as the principal directions
      auto conc = (*_concentration)[_qp]; 
      if (_concentration){
        //mooseWarning((*_concentration[_qp] >= 0.0, "Concentration less than 0");
        if (conc < 0.0)
                conc = 0.0;
    }

    Jg = 1.0 + _omega * (conc - _cref);
    ADReal Theta = _omega / Jg;
    // Volumetric growth strain rate 
      
      
    if (!_isotropic_swelling)
    {
        ADRankTwoTensor Fi;
        ADRealVectorValue m_v1;
        ADRealVectorValue m_v2;
        ADRealVectorValue m_v3;

        Fi = _Fe[_qp] * _Fp[_qp];    

        //Compute the first principal direction: normal direction
        if (_use_deformed_directions)
        {
          m_v1 = Fi.inverse() * _swell_normal[_qp];
        }else {
            m_v1 = _swell_normal[_qp];
        }
        m_v1 = m_v1 / m_v1.norm();

        //compute the second and third principal directions
        m_v3 = {0, 0, 1};
        m_v2 = m_v1.cross(m_v3);

        //compute the principal tensors
        ADRankTwoTensor S1, S2, S3;

        S1.vectorOuterProduct(m_v1, m_v1);
        S2.vectorOuterProduct(m_v2, m_v2);
        S3.vectorOuterProduct(m_v3, m_v3);

        ADRankTwoTensor Sg; 
        Sg = _alpha(0) * S1 + _alpha(1) * S2 + _alpha(2) * S3;

        _growth_tensor[_qp] = Sg;
        
        if (_rate_form )
        {
            ADReal epsilon_dot_g = Theta * (*_cdot)[_qp];
            ADRankTwoTensor Dgdt = epsilon_dot_g * Sg * _dt;
            // This computes the spectral decomposition of DgDt and computes
            // eigenvalues and eigenvectors... which in turn allows us to calculate 
            // the exponential of the DgDt tensor
            std::vector<ADReal> e_value(3);
            ADRankTwoTensor e_vector, N1, N2, N3;

            Dgdt.symmetricEigenvaluesEigenvectors(e_value, e_vector);

            N1.vectorOuterProduct(e_vector.column(0), e_vector.column(0));
            N2.vectorOuterProduct(e_vector.column(1), e_vector.column(1));
            N3.vectorOuterProduct(e_vector.column(2), e_vector.column(2));
            ADRankTwoTensor expDg = N1 * std::exp(e_value[0]) + 
                                  N2 * std::exp(e_value[1]) + 
                                  N3 * std::exp(e_value[2]);  
             Fg = expDg * _Fg_old[_qp];
        } else
        {
            ADReal lambda1, lambda2, lambda3;
            if (Jg > 0) 
            {
              if (_alpha(0) > 0) 
                  lambda1 = std::pow(Jg, _alpha(0));
              else
                  lambda1 = 1.0;

              if (_alpha(1) > 0) 
                  lambda2 = std::pow(Jg, _alpha(1));
              else
                  lambda2 = 1.0;

              if (_alpha(2) > 0) 
                  lambda3 = std::pow(Jg, _alpha(2));
              else
                  lambda3 = 1.0;
            } else {
                lambda1 = 1.0;
                lambda2 = 1.0;
                lambda3 = 1.0;
            }

            Fg  = lambda1 * S1 + lambda2 * S2 + lambda3 * S3; 
            
        }
        mooseAssert(Fg.det() > 0.0, "Det of Fg less than 0");
    } else
    {
        _growth_tensor[_qp] = _alpha(0) * _identity_two;
        Fg = std::pow(Jg, _alpha(0)) * _identity_two;
    }
       
}  

void
ADIsoTropicHyperViscoSwellingBase::computeNormaltoBoundaryandPrincipalDirections()
{
    
    if (_boundary_name.empty())
    {
        _swell_normal[_qp] = {0,1,0};
        return;
    }
// Define the minimal two distance from the _qp point to the two nodes at boundary
  Real minD1 = 1000;
  Real minD2 = 1000;
// The distance from the _qp point to the current node at boundary
  Real nrmQn = 0;
// The vectors from the _qp point to the two nodes at boundary that are shortest to the _qp point
  RealVectorValue vec1 = {0, 0, 0};
  RealVectorValue vec2 = {0, 0, 0};
// For temporary use
  RealVectorValue tmpVec;

// Get the node_id list at the interface boundary
  std::vector<dof_id_type> nodelist;
  for (auto & bnd_name : _boundary_name)
    nodelist = _mesh.getNodeList(_mesh.getBoundaryID(bnd_name));

// Compute all the distance, and find the two shortest-distance vector 
  Node * _node_ptr;
  for (auto i = 0; i < nodelist.size(); ++i)
  {
// Get the coordinate of nodal point at interface boundary
    _node_ptr = _mesh.getMesh().query_node_ptr(nodelist[i]);
    const Point & p = * _node_ptr;
    tmpVec   = {p(0), p(1), p(2)};
// Get the distance from the current qp point to the nodal point
    tmpVec  -= _q_point[_qp];
    nrmQn = tmpVec.norm();
// Find the shortest-distance vector
    if (nrmQn < minD1)
    {
// The 2nd_shortest takes the old 1st_shortest value
      minD2 = minD1;
      vec2  = vec1;
// The new 1st_shortest takes the current shortest value
      minD1 = nrmQn;
      vec1  = tmpVec;
    } 
// If the distance is larger than 2nd_shortest
    else if (nrmQn < minD2)
    {
      minD2 = nrmQn;
      vec2  = tmpVec;
    }

  }
 
// Compute the tangential vector connecting the two nodes
  // This current evaluation assumes that the 2 tangential directions are ortho-
  // normal. This is not strictly correct. So we need to compute prinicipal directions
  // based on normal to boundary in the current config and transform all vectors
  // to the intermediate space
  tmpVec = vec1 - vec2;
  tmpVec = tmpVec/tmpVec.norm();

// Compute the swell normals, make sure the direction is always point inward of Metal; only for 2D case
  tmpVec = {-tmpVec(1),tmpVec(0),tmpVec(2)};
  if (tmpVec(1) < 0)
    tmpVec = -tmpVec;
// Store the normals into material property
  _swell_normal[_qp] = tmpVec;
  
  // Debug printing of _qp principal directions for verification 
  // This will eventually be removed in favor of a more MOOSE based approach
#ifdef DEBUG
// Output result for debug usage
  std::ofstream myfile;
  myfile.open ("out.txt",std::ios::app);
  myfile << _q_point[_qp](0) << " " << _q_point[_qp](1) << " " << tmpVec(0) << " " << tmpVec(1) << std::endl;
  myfile.close();
#endif 

}
