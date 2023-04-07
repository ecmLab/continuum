
#include "SingleSE.h"
#include "Function.h"

registerADMooseObject("ecmElectrochemApp", SingleSE);

InputParameters
SingleSE::validParams()
{
   InputParameters params = ADMaterial::validParams();

// Add a required parameter.  If this isn't provided in the input file MOOSE will error.
  params.addRequiredParam<Real>("ionic_conductivity", "The ionic conductivity in bulk SE.");

// Add a required parameter.  If this isn't provided in the input file MOOSE will error.
  params.addRequiredParam<Real>("gb_conductivity", "The ionic conductivity in grain boundary.");

// Add a required parameter.  If this isn't provided in the input file MOOSE will error.
  params.addRequiredParam<Real>("metal_conductivity", "The electric conductivity in Li metal).");

// Add a required parameter.  If this isn't provided in the input file MOOSE will error.
    params.addRequiredParam<Real>("applied_current", "The applied current density from inlet boundary ($\\mathrm{A/cm^2}$).");

// Add a required parameter.  If this isn't provided in the input file MOOSE will error.
    params.addRequiredParam<Real>("exchange_current", "The exchange current density at the interface ($\\mathrm{A/cm^2}$).");

// Add a required parameter.  If this isn't provided in the input file MOOSE will error.
    params.addRequiredParam<Real>("reaction_rate", "The reaction_rate of Li+ to Li, unitless.");

// Add a required parameter.  If this isn't provided in the input file MOOSE will error.
    params.addRequiredParam<Real>("LiPotAnode", "The Li chemical potential in Anode, unit V.");

// Add a required parameter.  If this isn't provided in the input file MOOSE will error.
    params.addRequiredParam<Real>("LiPotCathode", "The Li chemical potential in Cathode, unit V.");
  
   return params;
}

SingleSE::SingleSE(const InputParameters & parameters)
  : ADMaterial(parameters),

    // Get the parameters from the input file
    _inIonicConductivity(getParam<Real>("ionic_conductivity")),
    _inGbConductivity(getParam<Real>("gb_conductivity")),
    _inMetalConductivity(getParam<Real>("metal_conductivity")),
    _inInletCurrent(getParam<Real>("applied_current")),
    _inExchangeCurrent(getParam<Real>("exchange_current")),
    _inReactionRate(getParam<Real>("reaction_rate")),
    _inLiPotAnode(getParam<Real>("LiPotAnode")),
    _inLiPotCathode(getParam<Real>("LiPotCathode")),

    // Declare material properties by getting a reference from the MOOSE Material system
    _ionic_conductivity(declareADProperty<Real>("ionic_conductivity")),
    _gb_conductivity(declareADProperty<Real>("gb_conductivity")),
    _metal_conductivity(declareADProperty<Real>("metal_conductivity")),
    _applied_current(declareADProperty<Real>("applied_current")),
    _exchange_current(declareADProperty<Real>("exchange_current")),
    _reaction_rate(declareADProperty<Real>("reaction_rate")),
    _LiPotAnode(declareADProperty<Real>("LiPotAnode")),
    _LiPotCathode(declareADProperty<Real>("LiPotCathode"))

{
}

void
SingleSE::computeQpProperties()
{
  _ionic_conductivity[_qp] = _inIonicConductivity;
  _gb_conductivity[_qp] = _inGbConductivity;
  _metal_conductivity[_qp] = _inMetalConductivity;
  _applied_current[_qp] = _inInletCurrent;
  _exchange_current[_qp] = _inExchangeCurrent;
  _reaction_rate[_qp] = _inReactionRate;
  _LiPotAnode[_qp] = _inLiPotAnode;
  _LiPotCathode[_qp] = _inLiPotCathode;
}
