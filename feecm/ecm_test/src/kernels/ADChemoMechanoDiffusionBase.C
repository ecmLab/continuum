

#include "ADChemoMechanoDiffusionBase.h"

//registerMooseObject("ecmApp", ADChemoMechanoDiffusion);
//registerMooseObject("ecmApp", ADChemoMechanoAnisoDiffusion);
//
//
//template<typename T>
//InputParameters
//ADChemoMechanoDiffusionTempl<T>::validParams()
//{
//  InputParameters params = ADMatDiffusionBase<T>::validParams();
//  params.addClassDescription("Diffusion equation kernel that takes an anisotropic diffusivity "
//                             "from a material property and "
//                             "computes the gradient of the chemical potential from "
//                             "elastic and growth contributions");
//  params.addCoupledVar("stress_based_chemical_potential",
//                          "Name of the variable for the stress_based_chemical_potential");
//
//  params.set<bool>("use_displaced_mesh") = false;
//  return params;
//}
//
//template <typename T>
//ADChemoMechanoDiffusionTempl<T>::ADChemoMechanoDiffusionTempl(const InputParameters & parameters)
//  : ADMatDiffusionBase<T>(parameters),
//    _mu_coupled(isCoupled("stress_based_chemical_potential")),
//    _mu_var(_mu_coupled ? coupled("stress_based_chemical_potential") : 0)
//{
//    if (_mu_coupled)
//        _grad_mu = &adCoupledGradient("stress_based_chemical_potential");
//
//}
//
//template<typename T>
//ADRealVectorValue
//ADChemoMechanoDiffusionTempl<T>::precomputeQpResidual()
//{
//    auto residual = ADMatDiffusionBase<T>::precomputeQpResidual();
//    if (_mu_coupled)
//        residual += _diffusivity[_qp] * _u[_qp] * (*_grad_mu)[_qp] ;
//
//    return residual;
//}
