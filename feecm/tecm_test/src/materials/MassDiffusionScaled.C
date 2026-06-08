// Copyright 2025, CEWLAB, All Rights Reserved
#include "MassDiffusionScaled.h"

registerMooseObject("tecm_testApp", MassDiffusion);

InputParameters
MassDiffusionScaled::validParams()
{
  InputParameters params = DualChemicalEnergyDensity::validParams();
  params.addRequiredParam<MaterialPropertyName>("mobility", "The mobility of the species");
  params.addParam<MaterialPropertyName>("Phi0", "The characteristic potential");
  params.addParam<MaterialPropertyName>("L0","The Characteristic Length");
  return params;
}

MassDiffusionScaled::MassDiffusionScaled(const InputParameters & parameters)
  : DualChemicalEnergyDensity(parameters),
  _M(getADMaterialProperty<Real>("mobility")),
  _phi0(getADMaterialProperty<Real>("Phi0")),
  _L0(getADMaterialProperty<Real>("L0"))
{
}

void
MassDiffusionScaled::computeQpProperties()
{
  _d_zeta_d_grad_mu[_qp] = -std::pow(_phi0[_qp],2)/std::pow(_L0[_qp],2)*_M[_qp] * _grad_mu[_qp];
  _zeta[_qp] = 0.5 * _d_zeta_d_grad_mu[_qp] * _grad_mu[_qp];
}
