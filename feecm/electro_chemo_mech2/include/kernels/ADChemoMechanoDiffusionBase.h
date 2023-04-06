    //* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
/* 
 * File:   ADChemoMechanoDiffusionBase.h
 * Author: srinath
 * 
 * Created on September 21, 2020, 10:57 AM
 */

#pragma once

#include "ADMatDiffusionBase.h"

template<typename T> 
class ADChemoMechanoDiffusionTempl : public ADMatDiffusionBase<T>
{
public: 
    static InputParameters validParams();

    ADChemoMechanoDiffusionTempl(const InputParameters & parameters);

protected:
    virtual ADRealVectorValue precomputeQpResidual() override;
    const bool _mu_coupled; 
    const unsigned int _mu_var;

    const ADVariableGradient * _grad_mu;
    using ADMatDiffusionBase<T> :: _qp;
    using ADMatDiffusionBase<T> :: _diffusivity;
    using ADMatDiffusionBase<T> :: _u;
    using ADMatDiffusionBase<T> :: adCoupledGradient;
    using ADMatDiffusionBase<T> :: coupled;
    using ADMatDiffusionBase<T> :: isCoupled;
    using ADMatDiffusionBase<T> :: getParam;
    Real _gas_constant; 
    Real _temperature; 

private:

};

template<typename T> 
InputParameters
ADChemoMechanoDiffusionTempl<T>::validParams()
{
  InputParameters params = ADMatDiffusionBase<T>::validParams();
  params.addClassDescription("Diffusion equation kernel that takes an anisotropic diffusivity "
                             "from a material property and "
                             "computes the gradient of the chemical potential from "
                             "elastic and growth contributions");
  params.addCoupledVar("stress_based_chemical_potential", 
                          "Name of the variable for the stress_based_chemical_potential");
  params.addParam<Real>("gas_constant", 8.314426, "Universal Gas Constant");
  params.addParam<Real>("temperature", 298, "temperature");

  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

template <typename T>
ADChemoMechanoDiffusionTempl<T>::ADChemoMechanoDiffusionTempl(const InputParameters & parameters)
  : ADMatDiffusionBase<T>(parameters),
    _mu_coupled(isCoupled("stress_based_chemical_potential")),
    _mu_var(_mu_coupled ? coupled("stress_based_chemical_potential") : 0), 
    _gas_constant(this-> template getParam<Real>("gas_constant")), 
    _temperature(this-> template getParam<Real>("temperature"))

{
    if (_mu_coupled)
        _grad_mu = &adCoupledGradient("stress_based_chemical_potential");

}

template<typename T> 
ADRealVectorValue
ADChemoMechanoDiffusionTempl<T>::precomputeQpResidual()
{
    auto residual = ADMatDiffusionBase<T>::precomputeQpResidual();
    if (_mu_coupled)
        residual += _diffusivity[_qp] * _u[_qp] * (*_grad_mu)[_qp] 
                / (_gas_constant * _temperature) ;

    return residual;
}