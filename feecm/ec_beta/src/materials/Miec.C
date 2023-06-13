
#include "Miec.h"
#include "Function.h"

registerADMooseObject("ecBetaApp", Miec);

InputParameters Miec::validParams()
{
  InputParameters params = ADMaterial::validParams();

// Add a required parameter.  If this isn't provided in the input file MOOSE will error.
  params.addRequiredParam<Real>("ionic_conductivity", "The ionic conductivity in SE.");

// Add a required parameter.  If this isn't provided in the input file MOOSE will error.
  params.addRequiredParam<Real>("metal_conductivity", "The electric conductivity in Li metal).");

// Add a required parameter.  If this isn't provided in the input file MOOSE will error.
  params.addRequiredParam<Real>("electronic_conductivity", "The electronic conductivity in SE ($\\mathrm{mS/cm}$).");

// Add a required parameter.  If this isn't provided in the input file MOOSE will error.
    params.addRequiredParam<Real>("inlet_current", "The inlet current density from inlet boundary ($\\mathrm{A/cm^2}$).");

// Add a required parameter.  If this isn't provided in the input file MOOSE will error.
    params.addRequiredParam<Real>("exchange_current", "The exchange current density at the interface ($\\mathrm{A/cm^2}$).");

// Add a required parameter.  If this isn't provided in the input file MOOSE will error.
    params.addRequiredParam<Real>("reaction_rate", "The reaction_rate of Li+ to Li, unitless.");

// Add a required parameter.  If this isn't provided in the input file MOOSE will error.
    params.addRequiredParam<Real>("electron_concentration", "The relative electron concentration in the SE, unitless.");

// Add a required parameter.  If this isn't provided in the input file MOOSE will error.
    params.addRequiredParam<Real>("LiPotAnode", "The Li chemical potential in Anode, unit V.");

// Add a required parameter.  If this isn't provided in the input file MOOSE will error.
    params.addRequiredParam<Real>("LiPotCathode", "The Li chemical potential in Cathode, unit V.");
 
  return params;
}

Miec::Miec(const InputParameters & parameters)
  : ADMaterial(parameters),

    // Get the parameters from the input file
    _inIonicConductivity(getParam<Real>("ionic_conductivity")),
    _inMetalConductivity(getParam<Real>("metal_conductivity")),
    _inElectronicConductivity(getParam<Real>("electronic_conductivity")),
    _inInletCurrent(getParam<Real>("inlet_current")),
    _inExchangeCurrent(getParam<Real>("exchange_current")),
    _inReactionRate(getParam<Real>("reaction_rate")),
    _inElectronConcentration(getParam<Real>("electron_concentration")),
    _inLiPotAnode(getParam<Real>("LiPotAnode")),
    _inLiPotCathode(getParam<Real>("LiPotCathode")),

    // Declare material properties by getting a reference from the MOOSE Material system
    _ionic_conductivity(declareADProperty<Real>("ionic_conductivity")),
    _metal_conductivity(declareADProperty<Real>("metal_conductivity")),
    _electronic_conductivity(declareADProperty<Real>("electronic_conductivity")),
    _inlet_current(declareADProperty<Real>("inlet_current")),
    _exchange_current(declareADProperty<Real>("exchange_current")),
    _reaction_rate(declareADProperty<Real>("reaction_rate")),
    _electron_concentration(declareADProperty<Real>("electron_concentration")),
    _LiPotAnode(declareADProperty<Real>("LiPotAnode")),
    _LiPotCathode(declareADProperty<Real>("LiPotCathode"))

{
}

void
Miec::computeQpProperties()
{
  _ionic_conductivity[_qp] = _inIonicConductivity;
  _metal_conductivity[_qp] = _inMetalConductivity;
  _electronic_conductivity[_qp] = _inElectronicConductivity;
  _inlet_current[_qp] = _inInletCurrent;
  _exchange_current[_qp] = _inExchangeCurrent;
  _reaction_rate[_qp] = _inReactionRate;
  _electron_concentration[_qp] = _inElectronConcentration;
  _LiPotAnode[_qp] = _inLiPotAnode;
  _LiPotCathode[_qp] = _inLiPotCathode;
}
