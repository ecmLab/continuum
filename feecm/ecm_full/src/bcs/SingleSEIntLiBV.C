
#include "SingleSEIntLiBV.h"

registerADMooseObject("electro_chemo_mechApp", SingleSEIntLiBV);

InputParameters
SingleSEIntLiBV::validParams()
{

  InputParameters params = ADIntegratedBC::validParams();
  params.addClassDescription("Compute the boundary condition of single conduction.");

// Add a coupled parameter: potLi
  params.addRequiredCoupledVar("potLi", "The potential of the other");

// Add a parameter with a default value; this value can be overridden in the input file.
  params.addParam<Real>(
        "F_RT",
        38.68,
        "The constant of F/RT, when T = 300K.");

  return params;
}

SingleSEIntLiBV::SingleSEIntLiBV(const InputParameters & parameters)
  : ADIntegratedBC(parameters),
  // Couple to the other potential
   _potLi(adCoupledValue("potLi")),

  // Get the parameters from the input file
   _exchange_current(getADMaterialProperty<Real>("exchange_current")),
   _reaction_rate(getADMaterialProperty<Real>("reaction_rate")),
   _F_RT(getParam<Real>("F_RT"))
{
}

ADReal
SingleSEIntLiBV::computeQpResidual()
{
  ADReal k1 = std::exp(_reaction_rate[_qp] * _F_RT * (_potLi[_qp] - _u[_qp]));
  ADReal k2 = std::exp(- (1 - _reaction_rate[_qp]) * _F_RT * (_potLi[_qp] - _u[_qp]));

  return _test[_i][_qp] * _exchange_current[_qp] * (k1 - k2);

}
