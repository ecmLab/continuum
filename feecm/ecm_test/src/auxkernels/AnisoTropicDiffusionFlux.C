/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   AnisoTropicDiffusionFlux.C
 * Author: srinath
 * 
 * Created on August 17, 2020, 4:57 PM
 */

#include "AnisoTropicDiffusionFlux.h"

registerMooseObject("electro_chemo_mechApp", AnisoTropicDiffusionFluxAux);

// defineLegacyParams(AnisoTropicDiffusionFluxAux);

InputParameters
AnisoTropicDiffusionFluxAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  MooseEnum component("x y z");
  params.addClassDescription("Compute components of flux vector for diffusion problems "
                             "with anisotropic diffusivity $(\\vv{J} = -D \\nabla C)$.");
  params.addRequiredParam<MooseEnum>("component", component, "The desired component of flux.");
  params.addRequiredCoupledVar("diffusion_variable", "The name of the variable");
  params.addParam<MaterialPropertyName>(
      "diffusivity",
      "diffusivity",
      "The name of the diffusivity material property that will be used in the flux computation.");

  return params;
}

AnisoTropicDiffusionFluxAux::AnisoTropicDiffusionFluxAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _component(getParam<MooseEnum>("component")),
    _grad_u(coupledGradient("diffusion_variable")),
    _diffusivity(getADMaterialProperty<RealTensorValue>("diffusivity"))
{
}

Real
AnisoTropicDiffusionFluxAux::computeValue()
{
  return MetaPhysicL::raw_value((-_diffusivity[_qp] * _grad_u[_qp])(_component));
}
