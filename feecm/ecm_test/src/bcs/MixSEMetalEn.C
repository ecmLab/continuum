
#include "MixSEMetalEn.h"

registerADMooseObject("ecmApp", MixSEMetalEn);

InputParameters
MixSEMetalEn::validParams()
{

  InputParameters params = ADIntegratedBC::validParams();
  params.addClassDescription("Compute the deposition boundary condition.");

// Add a coupled parameter: potMt
  params.addRequiredCoupledVar("potMt", "The potential of electron in metal of pore");

// Add a coupled parameter: potLi
  params.addRequiredCoupledVar("potLi", "The potential of Li+ in SE");

// Add a parameter with a default value; this value can be overridden in the input file.
  params.addParam<Real>(
        "F_RT",
        38.68,
        "The constant of F/RT, when T = 300K.");

  return params;
}

MixSEMetalEn::MixSEMetalEn(const InputParameters & parameters)
  : ADIntegratedBC(parameters),

  // Couple to the electron potential in metal of the pore
   _potMt(adCoupledValue("potMt")),

  // Get the gradient of the variable
   _potMt_gradient(adCoupledGradient("potMt")),

  // Couple to the Li potential
   _potLi(adCoupledValue("potLi")),

  // Get the parameters from the input file
   _metal_conductivity(getADMaterialProperty<Real>("metal_conductivity")),
   _exchange_current(getADMaterialProperty<Real>("exchange_current")),
   _reaction_rate(getADMaterialProperty<Real>("reaction_rate")),
   _electron_concentration(getADMaterialProperty<Real>("electron_concentration")),
   _F_RT(getParam<Real>("F_RT"))
{
}

ADReal
MixSEMetalEn::computeQpResidual()
{
  ADReal k1 = std::exp(_reaction_rate[_qp] * _F_RT * (_potLi[_qp] + _u[_qp] ));
  ADReal k2 = std::exp(- (1 - _reaction_rate[_qp]) * _F_RT * (_potLi[_qp] + _u[_qp] ));
  return _test[_i][_qp] * (-10000*_metal_conductivity[_qp] * _potMt_gradient[_qp] * _normals[_qp] + _electron_concentration[_qp] * _exchange_current[_qp] * (k1 - k2));
}
