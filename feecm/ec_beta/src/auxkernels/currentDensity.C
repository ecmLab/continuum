#include "currentDensity.h"

registerMooseObject("ecBetaApp", currentDensity);

InputParameters
currentDensity::validParams()
{
  InputParameters params = AuxKernel::validParams();

  // Declare the options for a MooseEnum.
  // These options will be presented to the user in Peacock and if something other than these
  // options is in the input file an error will be printed
  MooseEnum component("x y z");

  // Use the MooseEnum to add a parameter called "component"
  params.addRequiredParam<MooseEnum>("component", component, "The desired component of current density.");
  params.addRequiredParam<MaterialPropertyName>("conductivity", "The conductivity, in mS/cm.");

  // Add a "coupling paramater" to get a variable from the input file.
  params.addRequiredCoupledVar("potential", "The potential field.");

  return params;
}

currentDensity::currentDensity(const InputParameters & parameters)
  : AuxKernel(parameters),

    // Automatically convert the MooseEnum to an integer
    _component(getParam<MooseEnum>("component")),

    // Get the gradient of the variable
    _potential_gradient(coupledGradient("potential")),

    // Get the parameters from the input file
    _conductivity(parameters.get<MaterialPropertyName>("conductivity")),
    _conductivity_coef(getADMaterialProperty<Real>(_conductivity))

{
}

Real
currentDensity::computeValue()
{
  // Access the gradient of the potential at this quadrature point, then pull out the "component" of
  // it requested (x, y, z). Note, that getting a particular component of a gradient is done using
  // the parenthesis operator.

  Real k = 10.0 * MetaPhysicL::raw_value(_conductivity_coef[_qp]) * _potential_gradient[_qp](_component); // The prefactor 10 is due to unit conversion
  if(_conductivity == "ionic_conductivity") {
    k = -k;
  }
  return k;

}
