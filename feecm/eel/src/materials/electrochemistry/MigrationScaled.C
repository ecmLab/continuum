// Copyright 2023, UChicago Argonne, LLC All Rights Reserved
// License: L-GPL 3.0
#include "MigrationScaled.h"

registerMooseObject("EelApp", MigrationScaled);

InputParameters
MigrationScaled::validParams()
{
  InputParameters params = ElectroChemicalEnergyDensity::validParams();
  params.addClassDescription(
      params.getClassDescription() +
      " This class defines the electrochemical potential for the migration mechanism");
  params.addRequiredParam<MaterialPropertyName>("electric_conductivity",
                                                "The electric conductivity tensor");
  params.addRequiredParam<Real>("faraday_constant", "Faraday's constant");
  params.addParam<Real>("R", "Gas Constant");
  params.addParam<Real>("T", "Temperature");
  params.addParam<MaterialPropertyName>("L0", "Characteristic Length");
  params.addParam<MaterialPropertyName>("Phi0", "Characteristic potential");
  return params;
}

MigrationScaled::MigrationScaled(const InputParameters & parameters)
  : ElectroChemicalEnergyDensity(parameters),
    _sigma(getADMaterialProperty<Real>("electric_conductivity")),
    _F(getParam<Real>("faraday_constant")),
    _R(getParam<Real>("R")),
    _T(getParam<Real>("T")),
    _L0(getADMaterialProperty<Real>("L0")),
    _phi0(getADMaterialProperty<Real>("Phi0"))
{
}

void
MigrationScaled::computeQpProperties()
{
  _E[_qp] = _sigma[_qp] * _phi0[_qp]*_R*_T / (_F*_L0[_qp]*_L0[_qp]) * _grad_Phi[_qp] * _grad_mu[_qp];
  _d_E_d_grad_Phi[_qp] = _sigma[_qp]* _R*_T / (_F*_L0[_qp]) * _grad_mu[_qp];
  _d_E_d_grad_mu[_qp] = _sigma[_qp]*_phi0[_qp] / (_F*_L0[_qp]) * _grad_Phi[_qp];
}
