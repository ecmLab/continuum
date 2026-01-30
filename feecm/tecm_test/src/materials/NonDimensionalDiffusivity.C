// Copyright 2025, CEWLAB, All Rights Reserved
// License: L-GPL 3.0

#include "NonDimensionalDiffusivity.h"

registerMooseObject("tecm_testApp", NonDimensionalDiffusivity);

InputParameters
NonDimensionalDiffusivity::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription("Computes dimensionless diffusivity D̃ = D/D₀ using "
                            "characteristic scales from NonDimensionalParameters");

  params.addRequiredParam<Real>("D_dimensional", 
                               "Dimensional diffusivity [m²/s]");

  return params;
}

NonDimensionalDiffusivity::NonDimensionalDiffusivity(const InputParameters & parameters)
  : ADMaterial(parameters),
    _D_dimensional(getParam<Real>("D_dimensional")),
    _D0(getMaterialProperty<Real>("D0")),
    _diffusivity_nondim(declareADProperty<Real>("diffusivity_nondim")),
    _diffusivity_dim(declareADProperty<Real>("diffusivity_dim"))
{
  // Verify that NonDimensionalParameters material is present
  if (!hasMaterialProperty<Real>("D0"))
    paramError("D0", "NonDimensionalParameters material must be present to provide D0");
}

void
NonDimensionalDiffusivity::computeQpProperties()
{
  // Dimensional diffusivity (for reference/comparison)
  _diffusivity_dim[_qp] = _D_dimensional;
  
  // Dimensionless diffusivity: D̃ = D/D₀
  _diffusivity_nondim[_qp] = _D_dimensional / _D0[_qp];
}