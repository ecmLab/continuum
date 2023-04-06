/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   ADDiffusionFluxAux.C
 * Author: srinath
 *
 * Created on October 18, 2019, 8:50 AM
 */

#include "ADDiffusionFluxAux.h"

registerMooseObject("electro_chemo_mechApp", ADDiffusionFluxAux);

defineLegacyParams(ADDiffusionFluxAux);

InputParameters
ADDiffusionFluxAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  MooseEnum component("x y z");
  params.addClassDescription("Compute components of flux vector for diffusion problems "
                             "$(\\vv{J} = -D \\nabla C)$.");
  params.addRequiredParam<MooseEnum>("component", component, "The desired component of flux.");
  params.addRequiredCoupledVar("diffusion_variable", "The name of the variable");
  params.addRequiredParam<MaterialPropertyName>(
      "diffusivity",
      "The name of the diffusivity material property that will be used in the flux computation.");
  return params;
}

ADDiffusionFluxAux::ADDiffusionFluxAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _component(getParam<MooseEnum>("component")),
    _grad_u(coupledGradient("diffusion_variable")),
    _diffusion_coef(getADMaterialProperty<Real>("diffusivity"))
{
}

Real
ADDiffusionFluxAux::computeValue()
{
  return MetaPhysicL::raw_value(-_diffusion_coef[_qp]) * _grad_u[_qp](_component);
}
