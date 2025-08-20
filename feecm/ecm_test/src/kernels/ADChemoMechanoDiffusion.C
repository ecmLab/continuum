

#include "ADChemoMechanoDiffusion.h"


registerMooseObject("ecmApp", ADChemoMechanoDiffusion);

InputParameters
ADChemoMechanoDiffusion::validParams()
{
  auto params = ADChemoMechanoDiffusionTempl<Real>::validParams();
  params.addClassDescription(
      "Diffusion equation kernel that takes an isotropic diffusivity from a material property");

  return params;
}

ADChemoMechanoDiffusion::ADChemoMechanoDiffusion(const InputParameters & parameters)
  : ADChemoMechanoDiffusionTempl<Real>(parameters)
{
}
