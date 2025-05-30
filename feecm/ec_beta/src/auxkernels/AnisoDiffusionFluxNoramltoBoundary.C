/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   DiffusionFluxNormalToBoundaryAux.C
 * Author: srinath
 * 
 * Created on October 18, 2019, 8:50 AM
 */

#include "AnisoDiffusionFluxNoramltoBoundary.h"
#include "Assembly.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("ecBetaApp", AnisoDiffusionFluxNormalToBoundary);

// defineLegacyParams(AnisoDiffusionFluxNormalToBoundary);

InputParameters
AnisoDiffusionFluxNormalToBoundary::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Compute components of flux vector for diffusion problems "
                             "$(\\vv{J} = -D \\nabla C)$ normal to a boundary.");
  params.addRequiredCoupledVar("diffusion_variable", "The name of the variable");
  params.addRequiredParam<MaterialPropertyName>(
      "diffusivity",
      "The name of the diffusivity material property that will be used in the flux computation.");
  return params;
}

AnisoDiffusionFluxNormalToBoundary::AnisoDiffusionFluxNormalToBoundary(const InputParameters & parameters)
  : AuxKernel(parameters),
    _normals(_assembly.normals()),
    _grad_u(coupledGradient("diffusion_variable")),
    _diffusion_coef(getADMaterialProperty<RealTensorValue>("diffusivity"))
{
}

Real
AnisoDiffusionFluxNormalToBoundary::computeValue()
{
//    auto dim = 
//    RealVectorValue flux;
    auto flux = -_diffusion_coef[_qp] * _grad_u[_qp];
    return MetaPhysicL::raw_value(flux * _normals[_qp]);
}
