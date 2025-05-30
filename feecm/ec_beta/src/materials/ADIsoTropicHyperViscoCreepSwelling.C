//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADIsoTropicHyperViscoCreepSwelling.h"

#include "ElasticityTensorTools.h"

registerADMooseObject("ecBetaApp", ADIsoTropicHyperViscoSwellingCreep);


InputParameters
ADIsoTropicHyperViscoSwellingCreep::validParams()
{
	InputParameters params = ADIsoTropicHyperViscoSwellingBase::validParams();
    params.addClassDescription("Large strain Logarithmic strain based Istotropic Viscoplastic model "
                               "Based on Weber and Anand (1990) Computer Methods in App mech and engr"
                               "Hardening a function of strain rate");
    params.addRequiredParam<Real>("rate_exponent", "Rate Sensitivity parameter");
    params.addRequiredParam<Real>("initial_resistance", "Initial resistance in tension");
    params.addRequiredParam<Real>("hardening_modulus", "Hardening modulus in tension");
    params.addRequiredParam<Real>("saturation_resistance", "Saturation value of hardening in tension");
    params.addRequiredParam<Real>("hardening_exponent", "Hardening exponent (a)");
    params.addRequiredParam<Real>("pre_factor", "Pre factor for activation energy (A)");
    params.addRequiredParam<Real>("activation_energy", "Activation energy (Q)");
    params.addRequiredParam<Real>("gas_constant", "Universal gas constant (R)");
    params.addRequiredParam<Real>("temperature", "Temperature (T)");
    params.addRequiredParam<Real>("saturation_exponent", "Saturation Hardening exponent (n)");
        params.addDeprecatedParam<std::string>(
        "plastic_prepend",
        "",
        "String that is prepended to the plastic_strain Material Property",
        "This has been replaced by the 'base_name' parameter");
    params.set<std::string>("effective_inelastic_strain_name") = "effective_plastic_strain";
    return params;
}


ADIsoTropicHyperViscoSwellingCreep::ADIsoTropicHyperViscoSwellingCreep(const InputParameters & parameters)
  : ADIsoTropicHyperViscoSwellingBase(parameters),
    _mrate(getParam<Real>("rate_exponent")),
    _Y0(getParam<Real>("initial_resistance")),
    _H0(getParam<Real>("hardening_modulus")),
    _Ysat(getParam<Real>("saturation_resistance")),
    _ahard(getParam<Real>("hardening_exponent")),
    _A(getParam<Real>("pre_factor")),
    _Q(getParam<Real>("activation_energy")),
    _R(getParam<Real>("gas_constant")),
    _T(getParam<Real>("temperature")),
    _n(getParam<Real>("saturation_exponent")),
    _saturation_strength(declareADProperty<Real>(_base_name + _plastic_prepend + "saturation_strength")),
    _saturation_strength_old(getMaterialPropertyOld<Real>(_base_name + _plastic_prepend + "saturation_strength"))

{
    /* Convert to shear modulii like in Abaqus umat*/
     /* Convert to shear modulii like in Abaqus umat*/
    _prefactor = _A * std::exp(-_Q/(_R * _T));
    
    _epsilon_rate = _prefactor;
    _shear_rate = _epsilon_rate * std::sqrt(3.0);
    
    _shear_initial_resistance = _Y0 / std::sqrt(3.0);
    _shear_saturation = _Ysat / std::sqrt(3.0);
    _shear_initial_hardness = _H0 / 3.0;
    _check_range = true;   

}

void
ADIsoTropicHyperViscoSwellingCreep::initQpStatefulProperties()
{
   _hardening_variable[_qp] = _shear_initial_hardness;
  _strength_variable[_qp] = _shear_initial_resistance;
  _saturation_strength[_qp] = _shear_saturation;
  _yield_strength[_qp] = _Y0;
  ADIsoTropicHyperViscoSwellingBase::initQpStatefulProperties();
}

void
ADIsoTropicHyperViscoSwellingCreep::propagateQpStatefulProperties()
{
  _hardening_variable[_qp] = _hardening_variable_old[_qp];
  _strength_variable[_qp] = _strength_variable_old[_qp];
  _saturation_strength[_qp] = _saturation_strength_old[_qp];

  ADIsoTropicHyperViscoSwellingBase::propagateQpStatefulProperties();
}

void
ADIsoTropicHyperViscoSwellingCreep::computeStressInitialize(
    const ADReal & effective_trial_stress, const ADRankFourTensor & elasticity_tensor)
{
  _yield_strength[_qp] = _strength_variable[_qp] * std::sqrt(3.0);
}

void
ADIsoTropicHyperViscoSwellingCreep::iterationFinalize(ADReal scalar)
{
    if (scalar >= 0.0)
        _strength_variable[_qp] = computeHardeningValue(scalar);
        
    else
        _strength_variable[_qp] = _strength_variable_old[_qp];
}

ADReal
ADIsoTropicHyperViscoSwellingCreep::computeResidual(const ADReal & effective_trial_stress, const ADReal & scalar)
{
  ADReal residual = 0.0;
    
  if (scalar < 0) 
      return  effective_trial_stress/_three_shear_modulus/3.0;
      
  _strength_variable[_qp] = computeHardeningValue(scalar);
  residual = effective_trial_stress - _dt * scalar *_three_shear_modulus/3.0
            - _strength_variable[_qp] * (std::pow((scalar/_shear_rate), _mrate));
  residual /= _three_shear_modulus/3.0;
  return residual;
}

ADReal
ADIsoTropicHyperViscoSwellingCreep::computeDerivative(const ADReal & /*effective_trial_stress*/, const ADReal & scalar)
{
    ADReal derivative = 1.0;

    if (scalar < 1e-8)
            return -_dt;
    if (scalar > 0) {        
        auto hslope = computeHardeningDerivative(scalar);
        auto fac1 = hslope * std::pow((scalar/_shear_rate), _mrate);
        auto fac2 = _strength_variable[_qp] * (_mrate/_shear_rate) * 
                   std::pow((scalar / _shear_rate), _mrate-1.0);
        derivative = -_three_shear_modulus/3.0 * _dt - fac1 -fac2;
    } else
        derivative = -_three_shear_modulus/3.0 * _dt;
    derivative /= _three_shear_modulus/3.0;
    return derivative;
}

void
ADIsoTropicHyperViscoSwellingCreep::computeStressFinalize(const ADRankTwoTensor & plastic_strain_increment)
{
  _plastic_strain[_qp] += plastic_strain_increment;
  _yield_strength[_qp] = _strength_variable[_qp]*std::sqrt(3.0);
}

ADReal
ADIsoTropicHyperViscoSwellingCreep::computeHardeningValue(const ADReal & scalar)
{
    ADReal ssat;
    if (scalar > 0)
        ssat = _shear_saturation * std::pow(scalar/_shear_rate , _n);
    else
        ssat = _shear_saturation; 
    _saturation_strength[_qp] = ssat;
    
    auto temp = 1.0 - _strength_variable_old[_qp] / ssat; 
    auto sg = (temp < 0 ? -1.0 : 1.0);
    _hardening_variable[_qp] = _shear_initial_hardness 
                        * std::pow(std::abs(temp), _ahard) * sg;
    return _strength_variable_old[_qp] + _dt * _hardening_variable[_qp] * scalar;
}

ADReal
ADIsoTropicHyperViscoSwellingCreep::computeHardeningDerivative(
    const ADReal & /*scalar*/)
{
 return _hardening_variable[_qp] * _dt;
}

Real
ADIsoTropicHyperViscoSwellingCreep:: computeTimeStepLimit()
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
