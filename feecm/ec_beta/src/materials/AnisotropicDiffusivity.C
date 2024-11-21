#include "AnisotropicDiffusivity.h"
#include "Function.h"

// Register the material with MOOSE application "ecBetaApp"
registerADMooseObject("ecBetaApp", AnisotropicDiffusivity);

// Define valid parameters for AnisotropicDiffusivity material
InputParameters
AnisotropicDiffusivity::validParams()

{
    InputParameters params = ADMaterial::validParams();
    params.addRequiredParam<Real>("k0", "Base diffusivity");
    params.addRequiredParam<Real>("w", "Amplitude of the cosine term");
    params.addRequiredParam<Real>("lambda", "Wavelength parameter");
    params.addRequiredCoupledVar("gradient_variable", "The variable from which to compute the gradient component");

    return params;
}

// Constructor implementation
AnisotropicDiffusivity::AnisotropicDiffusivity(const InputParameters & parameters)
  : ADMaterial(parameters),
    _k0(getParam<Real>("k0")),
    _w(getParam<Real>("w")),
    _lambda(getParam<Real>("lambda")),
    _gradient(coupledGradient("gradient_variable")),
    _k(declareADProperty<Real>("diffusivity"))
{
}
void AnisotropicDiffusivity::computeQpProperties()
{
    // Get the gradient components in x and y directions
    const Real dop_dx = _gradient[_qp](0); // Gradient in x-direction
    const Real dop_dy = _gradient[_qp](1); // Gradient in y-direction

    // Calculate theta using atan2 to handle quadrant issues correctly
    Real theta = atan2(dop_dy, dop_dx);

    _k[_qp] = _k0 * pow((1 + _w * cos(_lambda * theta)), 2);
}
