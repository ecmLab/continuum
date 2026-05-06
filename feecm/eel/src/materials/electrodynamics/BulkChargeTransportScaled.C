// Copyright 2023, UChicago Argonne, LLC All Rights Reserved
// License: L-GPL 3.0
#include "BulkChargeTransportScaled.h"

registerMooseObject("EelApp", BulkChargeTransportScaled);

InputParameters
BulkChargeTransportScaled::validParams()
{
  InputParameters params = ElectricalEnergyDensity::validParams();
  params.addClassDescription(
      params.getClassDescription() +
      " This class defines the electrical potential for charge transfer in the bulk");
  params.addRequiredParam<MaterialPropertyName>("electric_conductivity",
                                                "The electric conductivity tensor");
  params.addParam<MaterialPropertyName>("Phi0", "The characteristic potential");
  params.addParam<MaterialPropertyName>("L0","The Characteristic Length");
  params.addParam<VariableName>("temperature", "The temperature");
  return params;
}

BulkChargeTransportScaled::BulkChargeTransportScaled(const InputParameters & parameters)
  : ElectricalEnergyDensity(parameters),
    _sigma(getADMaterialProperty<Real>("electric_conductivity")),
    _phi0(getADMaterialProperty<Real>("Phi0")),
    _L0(getADMaterialProperty<Real>("L0")),
    _d_E_d_lnT(isParamValid("temperature")
                   ? &declarePropertyDerivative<Real, true>(
                         _energy_name, "ln(" + getParam<VariableName>("temperature") + ")")
                   : nullptr)
{
}

void
BulkChargeTransportScaled::computeQpProperties()
{
  _d_E_d_grad_Phi[_qp] = std::pow(_phi0[_qp],2)/std::pow(_L0[_qp],2)*_sigma[_qp] * _grad_Phi[_qp];
  _E[_qp] = 0.5 * _d_E_d_grad_Phi[_qp] * _grad_Phi[_qp];

  if (_d_E_d_lnT)
    (*_d_E_d_lnT)[_qp] = _d_E_d_grad_Phi[_qp] * _grad_Phi[_qp];
}
