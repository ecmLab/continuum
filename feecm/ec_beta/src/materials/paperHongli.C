#include "paperHongli.h"
#include "Function.h"

registerMooseObject("ecBetaApp", paperHongli);

InputParameters
paperHongli::validParams()
{
  InputParameters params = Material::validParams();

  // Add a required parameter.  If this isn't provided in the input file MOOSE will error.
  params.addRequiredCoupledVar("cLi", "The Li content in the graphite.");
  params.addRequiredParam<Real>("c1", "The first constant in the diffusivity-cLi relation");
  params.addRequiredParam<Real>("c2", "The second constant in the diffusivity-cLi relation");

  return params;
}

paperHongli::paperHongli(const InputParameters & parameters)
  : Material(parameters),
    
    // Couple to the Li concentration
    _cLi(adCoupledValue("cLi")),
    // Get parameters from the input file
    _c1(getParam<Real>("c1")),
    _c2(getParam<Real>("c2")),

    // Declare material properties by getting a reference from the MOOSE Material system
   _diffusivity(declareADProperty<Real>("diffusivity")),
   _ionic_conductivity(declareADProperty<Real>("ionic_conductivity"))
{
}

void
paperHongli::computeQpProperties()
{
//     	_diffusivity[_qp] = _c1 * std::exp(_c2*6*_cLi[_qp]);
     	_diffusivity[_qp] = _c1 * std::pow(10, _c2*6*_cLi[_qp]);
        _ionic_conductivity[_qp] = 0.11*_diffusivity[_qp];
}
