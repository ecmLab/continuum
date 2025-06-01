

#include "ADChemoMechanoAnisoDiffusion.h"


registerMooseObject("ecmApp", ADChemoMechanoAnsioDiffusion);

InputParameters
ADChemoMechanoAnsioDiffusion::validParams()
{
  auto params = ADChemoMechanoDiffusionTempl<RealTensorValue>::validParams();
  params.addClassDescription(
      "Diffusion equation kernel that takes an isotropic diffusivity from a material property");
  return params;
}

ADChemoMechanoAnsioDiffusion::ADChemoMechanoAnsioDiffusion(const InputParameters & parameters)
  : ADChemoMechanoDiffusionTempl<RealTensorValue>(parameters)

{
}
