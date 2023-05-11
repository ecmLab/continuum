
#include "ionicSE.h"
#include "Function.h"

registerADMooseObject("ecmElectrochemApp", ionicSE);

InputParameters
ionicSE::validParams()
{
   InputParameters params = ADMaterial::validParams();

////**** Bulk properties of the SE
// Add bulk ionic conductivity as a required parameter.  If this isn't provided in the input file MOOSE will error.
  params.addRequiredParam<Real>("ionic_conductivity", "The ionic conductivity in bulk SE.");
// Add grain boundary conductivity with a default value; this value can be overridden in the input file.
  params.addParam<Real>("gb_conductivity",0.0, "The ionic conductivity in grain boundary.");

////**** Interface properties of the SE with electrode
// Add the exchange current as a required parameter.  If this isn't provided in the input file MOOSE will error.
//    params.addRequiredParam<Real>("exchange_current", "The exchange current density at the interface in mA/cm^2.");
// Add reaction rate with a default value; this value can be overridden in the input file.
//    params.addParam<Real>("reaction_rate",0.5, "The reaction_rate of Li+ to Li, unitless.");
  
   return params;
}

ionicSE::ionicSE(const InputParameters & parameters)
  : ADMaterial(parameters),

    // Get the parameters from the input file
    _inIonicConductivity(getParam<Real>("ionic_conductivity")),
    _inGbConductivity(getParam<Real>("gb_conductivity")),
//    _inExchangeCurrent(getParam<Real>("exchange_current")),
//    _inReactionRate(getParam<Real>("reaction_rate")),

    // Declare material properties by getting a reference from the MOOSE Material system
    _ionic_conductivity(declareADProperty<Real>("ionic_conductivity")),
    _gb_conductivity(declareADProperty<Real>("gb_conductivity"))
//    _exchange_current(declareADProperty<Real>("exchange_current")),
//    _reaction_rate(declareADProperty<Real>("reaction_rate"))

{
}

void
ionicSE::computeQpProperties()
{
  _ionic_conductivity[_qp] = _inIonicConductivity;
  _gb_conductivity[_qp] = _inGbConductivity;
//  _exchange_current[_qp] = _inExchangeCurrent;
//  _reaction_rate[_qp] = _inReactionRate;
}
