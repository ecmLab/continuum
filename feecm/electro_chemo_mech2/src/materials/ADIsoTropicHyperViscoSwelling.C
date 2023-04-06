//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADIsoTropicHyperViscoSwelling.h"

#include "ElasticityTensorTools.h"

registerADMooseObject("electro_chemo_mechApp", ADIsoTropicHyperViscoSwelling);

InputParameters ADIsoTropicHyperViscoSwelling::validParams()
{
	InputParameters params = ADIsoTropicHyperViscoSwellingBase::validParams();
    params.addClassDescription("Large strain Logarithmic strain based Istotropic Viscoplastic model "
                               "Based on Weber and Anand (1990) Computer Methods in App mech and engr"
                               "Hardening a function of strain rate");
    params.addRequiredParam<Real>("rate_exponent", "Rate Sensitivity parameter");
    params.addRequiredParam<Real>("reference_strain_rate", "Reference strain rate");
    params.addRequiredParam<Real>("initial_resistance", "Initial resistance in tension");
    params.addRequiredParam<Real>("hardening_modulus", "Hardening modulus in tension");
    params.addRequiredParam<Real>("saturation_resistance", "Saturation value of hardening in tension");
    params.addRequiredParam<Real>("hardening_exponent", "Hardening exponent");
        params.addDeprecatedParam<std::string>(
        "plastic_prepend",
        "",
        "String that is prepended to the plastic_strain Material Property",
        "This has been replaced by the 'base_name' parameter");
    params.set<std::string>("effective_inelastic_strain_name") = "effective_plastic_strain";

    return params;
}

ADIsoTropicHyperViscoSwelling::ADIsoTropicHyperViscoSwelling(const InputParameters & parameters)
  : ADIsoTropicHyperViscoSwellingBase(parameters),
    _mrate(getParam<Real>("rate_exponent")),
    _epsilon_rate(getParam<Real>("reference_strain_rate")),
    _Y0(getParam<Real>("initial_resistance")),
    _H0(getParam<Real>("hardening_modulus")),
    _Ysat(getParam<Real>("saturation_resistance")),
    _ahard(getParam<Real>("hardening_exponent"))
{
    /* Convert to shear modulii like in Abaqus umat*/
     /* Convert to shear modulii like in Abaqus umat*/
    _shear_rate = _epsilon_rate * std::sqrt(3.0);
    _shear_initial_resistance = _Y0 / std::sqrt(3.0);
    _shear_saturation = _Ysat / std::sqrt(3.0);
    _shear_initial_hardness = _H0 / 3.0;
    _check_range = true;   
//    _line_search = false;

    /*--------------------------------------------*/
//    _shear_rate = _epsilon_rate;
//    _shear_initial_resistance = _Y0;
//    _shear_saturation = _Ysat;
//    _shear_initial_hardness = _H0;

}



void
ADIsoTropicHyperViscoSwelling::initQpStatefulProperties()
{
   _hardening_variable[_qp] = _shear_initial_hardness;
  _strength_variable[_qp] = _shear_initial_resistance;
  _yield_strength[_qp] = _Y0;
  ADIsoTropicHyperViscoSwellingBase::initQpStatefulProperties();
}


void
ADIsoTropicHyperViscoSwelling::propagateQpStatefulProperties()
{
  _hardening_variable[_qp] = _hardening_variable_old[_qp];
  _strength_variable[_qp] = _strength_variable_old[_qp];
  ADIsoTropicHyperViscoSwellingBase::propagateQpStatefulProperties();
}


void
ADIsoTropicHyperViscoSwelling::computeStressInitialize(const ADReal & effective_trial_stress, const ADRankFourTensor & elasticity_tensor)
{
  
  _hardening_variable[_qp] = _shear_initial_hardness 
                            * std::pow((1.0 - (_strength_variable_old[_qp] / _shear_saturation)), _ahard);
  _yield_strength[_qp] = _strength_variable[_qp] * std::sqrt(3.0);
}



void
ADIsoTropicHyperViscoSwelling::iterationFinalize(ADReal scalar)
{
    if (scalar >= 0.0)
        _strength_variable[_qp] = _strength_variable_old[_qp] + _dt * _hardening_variable[_qp] * scalar;
        
    else
        _strength_variable[_qp] = _strength_variable_old[_qp];
}


ADReal
ADIsoTropicHyperViscoSwelling::computeResidual(const ADReal & effective_trial_stress, const ADReal & scalar)
{
  ADReal residual = 0.0;

//    if (scalar < 1.0e-9)
//            return 0.0;
    
    if (scalar < 0) 
        return  effective_trial_stress/_three_shear_modulus/3.0;

    _strength_variable[_qp] = computeHardeningValue(scalar);
    residual = effective_trial_stress - _dt * scalar *_three_shear_modulus/3.0
              - _strength_variable[_qp] * (std::pow((scalar/_shear_rate), _mrate));
    residual /= _three_shear_modulus/3.0;
  return residual;
}


ADReal
ADIsoTropicHyperViscoSwelling::computeDerivative(const ADReal & /*effective_trial_stress*/, const ADReal & scalar)
{
    ADReal derivative = 1.0;
//    if (scalar < 1e-9)
//        return 1.0;
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
ADIsoTropicHyperViscoSwelling::computeStressFinalize(const ADRankTwoTensor & plastic_strain_increment)
{
  _plastic_strain[_qp] += plastic_strain_increment;
  _yield_strength[_qp] = _strength_variable[_qp]*std::sqrt(3.0);
}


ADReal
ADIsoTropicHyperViscoSwelling::computeHardeningValue(const ADReal & scalar)
{
   return _strength_variable_old[_qp] + _dt * _hardening_variable[_qp] * scalar;

}

//template <>
//DualReal
//ADIsoTropicHyperViscoSwelling<JACOBIAN>::computeHardeningValue(const DualReal & scalar)
//{
//    return _strength_variable_old[_qp] + _dt * _hardening_variable[_qp] * scalar;
//}


ADReal
ADIsoTropicHyperViscoSwelling::computeHardeningDerivative(const ADReal & /*scalar*/)
{
 return _hardening_variable[_qp] * _dt;
}




Real
ADIsoTropicHyperViscoSwelling:: computeTimeStepLimit()
{
    Real scalar_inelastic_strain_incr;
    
    scalar_inelastic_strain_incr = std::abs (MetaPhysicL::raw_value(_effective_inelastic_strain[_qp] -
            _effective_inelastic_strain_old[_qp]))/ _max_inelastic_increment;
//    if (MooseUtils::absoluteFuzzyEqual(scalar_inelastic_strain_incr, 0.0))
//        return std::numeric_limits<Real>::max();

//    return _dt * _max_inelastic_increment / scalar_inelastic_strain_incr;
    if (scalar_inelastic_strain_incr <= 0.5)
        return _dt * 1.5;
    else if (scalar_inelastic_strain_incr > 0.5 && scalar_inelastic_strain_incr <= 0.8)
        return _dt * 1.25;
    else if (scalar_inelastic_strain_incr > 0.8 && scalar_inelastic_strain_incr <= 1.25)
        return _dt * 0.75;
    else
        return _dt * 0.5;
}
