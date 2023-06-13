#include "ChargedTransport.h"

registerADMooseObject("ecBetaApp", ChargedTransport);

InputParameters
ChargedTransport::validParams()
{
    InputParameters params = ADKernel::validParams();
    params.addClassDescription("Compute the laplacian equation.");
    params.addRequiredParam<MaterialPropertyName>("diffusivity", "The diffusivity coefficient.");
    return params;
}

ChargedTransport::ChargedTransport(const InputParameters & parameters)
  : ADKernel(parameters),

// Get the parameters from the input file
    _diffusivity(parameters.get<MaterialPropertyName>("diffusivity")),
    _diffusivity_coef(getADMaterialProperty<Real>(_diffusivity))
{
}

ADReal
ChargedTransport::computeQpResidual()
{
  ADReal k = 10.0 * MetaPhysicL::raw_value(_diffusivity_coef[_qp]) * _grad_test[_i][_qp] * _grad_u[_qp]; // The prefactor 10.0 is due to unit conversion

  if(_diffusivity == "electronic_conductivity") {
    k = -k;
  }
  return k;
}

