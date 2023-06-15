
#include "Miec.h"
#include "Function.h"

registerADMooseObject("ecBetaApp", Miec);

InputParameters
Miec::validParams()
{
   InputParameters params = ADMaterial::validParams();

////**** Bulk properties of the Miec
// Add bulk ionic conductivity as a required parameter.  If this isn't provided in the input file MOOSE will error.
  params.addRequiredParam<Real>("ionic_conductivity", "The ionic conductivity in bulk miec.");
// Add bulk electronic conductivity as a required parameter.  If this isn't provided in the input file MOOSE will error.
  params.addRequiredParam<Real>("electronic_conductivity", "The ionic conductivity in bulk miec.");
// Add grain boundary conductivity with a default value; this value can be overridden in the input file.
  params.addParam<Real>("gb_conductivity",0.0, "The ionic conductivity in grain boundary.");

   return params;
}

Miec::Miec(const InputParameters & parameters)
  : ADMaterial(parameters),

    // Get the parameters from the input file
    _inIonicConductivity(getParam<Real>("ionic_conductivity")),
    _inElectronicConductivity(getParam<Real>("electronic_conductivity")),
    _inGbConductivity(getParam<Real>("gb_conductivity")),

    // Declare material properties by getting a reference from the MOOSE Material system
    _ionic_conductivity(declareADProperty<Real>("ionic_conductivity")),
    _electronic_conductivity(declareADProperty<Real>("electronic_conductivity")),
    _gb_conductivity(declareADProperty<Real>("gb_conductivity"))
{
}

void
Miec::computeQpProperties()
{
  _ionic_conductivity[_qp] = _inIonicConductivity;
  _electronic_conductivity[_qp] = _inElectronicConductivity;
  _gb_conductivity[_qp] = _inGbConductivity;
}
