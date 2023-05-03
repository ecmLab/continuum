
#include "NernstPlanckConvection.h"

registerADMooseObject("ecmElectrochemApp", NernstPlanckConvection);

InputParameters
NernstPlanckConvection::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addRequiredCoupledVar("Voltage", "The variable representing the voltage.");
  params.addRequiredParam<MaterialPropertyName>("diffusivity", "The diffusivity coefficient.");
  params.addParam<Real>("zIons", 1, "The charge state of the ions, default be positive 1");

  //params.addRequiredParam<Real>("Q_asterik", "The heat of transport of Cu (Q*) in Sn at T in J/mol");
  // Add a parameter with a default value.  This value can be overriden in the input file.
  params.addParam<Real>("F_RT", 38.68, "The constant of F/RT,in unit 1/V, when T = 300K.");

  return params;
}

NernstPlanckConvection::NernstPlanckConvection(const InputParameters & parameters)
   : Kernel(parameters),
    // Save off the coupled variable identifier for use in
    // computeQpOffDiagJacobian
    _V_var(coupled("Voltage")),
    // Save off the coupled value for use in Residual
    _V(coupledValue("Voltage")),
    // Couple to the gradient of the pressure
    _grad_V(coupledGradient("Voltage")),
    // Grab necessary material properties
    _diffusivity(getMaterialProperty<Real>("diffusivity")),
    _zIons(getParam<Real>("zIons")),
    //_Qh(getParam<Real>("Q_asterik")),
    _F_RT(getParam<Real>("F_RT"))
{
}

Real
NernstPlanckConvection::computeQpResidual()
{
  RealVectorValue ion_velocity = _F_RT * _diffusivity[_qp] * _grad_V[_qp];
  return -_test[_i][_qp] * _zIons * ion_velocity * _grad_u[_qp];
}

Real
NernstPlanckConvection::computeQpJacobian()
{
  RealVectorValue ion_velocity = _F_RT * _diffusivity[_qp] * _grad_V[_qp];
  return -_test[_i][_qp] * _zIons * ion_velocity * _grad_phi[_j][_qp];
}

Real
NernstPlanckConvection::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _V_var)
  {
   return  -_test[_i][_qp] * _zIons * _F_RT * _diffusivity[_qp] * _grad_phi[_j][_qp] * _grad_u[_qp];
  }

  return 0.0;
}
