
#include "NernstPlanckConvectionEtime.h"

registerADMooseObject("ecBetaApp", NernstPlanckConvectionEtime);

InputParameters
NernstPlanckConvectionEtime::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addRequiredCoupledVar("Efield", "The variable representing the voltage.");
  params.addRequiredParam<MaterialPropertyName>("diffusivity", "The diffusivity coefficient.");
  params.addParam<Real>("zIons", 1, "The charge state of the ions, default be positive 1");
  params.addParam<Real>("F_RT", 38.68, "The constant of F/RT,in unit 1/V, when T = 300K.");
  params.addParam<Real>("scale", 1, "The constant of F/RT,in unit 1/V, when T = 300K.");

  return params;
}

NernstPlanckConvectionEtime::NernstPlanckConvectionEtime(const InputParameters & parameters)
   : Kernel(parameters),
    _Efield(coupledValue("Efield")),
    _diffusivity(getMaterialProperty<Real>("diffusivity")),
    _zIons(getParam<Real>("zIons")),
    _F_RT(getParam<Real>("F_RT")),
    _scale(getParam<Real>("scale"))
{
}

Real
NernstPlanckConvectionEtime::computeQpResidual()
{
  return -_scale*_test[_i][_qp] * _zIons * _F_RT * _diffusivity[_qp] * _u[_qp]*_Efield[_qp];
}

Real
NernstPlanckConvectionEtime::computeQpJacobian()
{
  return -_scale*_test[_i][_qp] * _zIons * _F_RT * _diffusivity[_qp] * _phi[_j][_qp]*_Efield[_qp];
}

Real
NernstPlanckConvectionEtime::computeQpOffDiagJacobian(unsigned int jvar)
{
  return 0.0;
}
