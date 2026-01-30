// Copyright 2025, CEWLAB, All Rights Reserved
// License: L-GPL 3.0

#include "NonDimensionalDiffusion.h"

registerMooseObject("tecm_testApp", NonDimensionalDiffusion);

InputParameters
NonDimensionalDiffusion::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Non-dimensional diffusion kernel for the TECM framework. "
                            "Implements -∇̃ · (D̃ ∇̃ũ) where all quantities are dimensionless.");

  params.addRequiredParam<MaterialPropertyName>(
      "diffusivity_nondim", 
      "Dimensionless diffusivity D̃ = D/D₀");

  params.addParam<bool>("use_nondimensional_scaling", 
                       true, 
                       "Whether to use non-dimensional scaling (true) or dimensional (false)");
                       
  params.addParam<MaterialPropertyName>(
      "diffusivity_dim",
      "Optional dimensional diffusivity for comparison mode");

  return params;
}

NonDimensionalDiffusion::NonDimensionalDiffusion(const InputParameters & parameters)
  : ADKernel(parameters),
    _diffusivity_nondim(getADMaterialProperty<Real>("diffusivity_nondim")),
    _L0(getMaterialProperty<Real>("L0")),
    _use_nondimensional_scaling(getParam<bool>("use_nondimensional_scaling")),
    _diffusivity_dim(isParamValid("diffusivity_dim") 
                     ? &getADMaterialProperty<Real>("diffusivity_dim") 
                     : nullptr)
{
  // Verify that NonDimensionalParameters material is present
  if (!hasMaterialProperty<Real>("L0"))
    paramError("L0", "NonDimensionalParameters material must be present to provide characteristic scales");
}

ADReal
NonDimensionalDiffusion::computeQpResidual()
{
  if (_use_nondimensional_scaling)
  {
    // Non-dimensional formulation: (D̃ ∇̃ũ, ∇̃φ)
    // The gradient operators are already dimensionless in the weak form
    // since we're solving for dimensionless variables
    return _diffusivity_nondim[_qp] * _grad_test[_i][_qp] * _grad_u[_qp];
  }
  else
  {
    // Dimensional formulation for comparison: (D ∇u, ∇φ) 
    if (_diffusivity_dim)
      return (*_diffusivity_dim)[_qp] * _grad_test[_i][_qp] * _grad_u[_qp];
    else
      mooseError("diffusivity_dim must be provided when use_nondimensional_scaling = false");
  }
}