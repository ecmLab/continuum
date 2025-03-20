#include "ADMigrationDiffusion.h"

registerADMooseObject("ecBetaApp", ADMigrationDiffusion);

InputParameters
ADMigrationDiffusion::validParams()
{
    InputParameters params = ADKernel::validParams();
    params.addClassDescription("Compute the the transport of ions due to migration and diffusion");
    params.addRequiredCoupledVar("c", "The Concentration field");
    params.addRequiredParam<MaterialPropertyName>("conductivity", "The conductivity coefficient.");
    params.addParam<Real>("c0", "The reference concentration");
    params.addParam<Real>("zIons", 1, "The charge state of the ions, default be positive 1");
    params.addParam<Real>("F_RT", 36.41, "The Value of F/RT at T=300K");
    params.addParam<Real>("scale", 1, "The constant of F/RT,in unit 1/V, when T = 300K.");
    return params;
}

ADMigrationDiffusion::ADMigrationDiffusion(const InputParameters & parameters)
  : ADKernel(parameters),
  _c(adCoupledValue("c")),
  _grad_c(adCoupledGradient("c")),
  _conductivity(getMaterialProperty<Real>("conductivity")),
  _c0(getParam<Real>("c0")),
  _zIons(getParam<Real>("zIons")),
  _F_RT(getParam<Real>("F_RT")),
  _scale(getParam<Real>("scale"))
{
}

ADReal
ADMigrationDiffusion::computeQpResidual()
{


  ADRealVectorValue term1 = _c0 / _c[_qp] * _grad_c[_qp];
  ADRealVectorValue term2 = _F_RT * _grad_u[_qp];
  return _scale*_conductivity[_qp]/_F_RT * _grad_test[_i][_qp] * (term1 + _zIons*term2);

}
