// Created by: Zeeshan Ahmad

#include "ChemPotValueHC.h"

registerMooseObject("edl_solidApp", ChemPotValueHC);

InputParameters
ChemPotValueHC::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Kernel to set the value of mu equal to a target value.");
  params.addRequiredParam<Real>("mutarget", "The target value of the chemical potential");
  params.addRequiredParam<MaterialPropertyName>("mu0", "The name of the first chemical potential");
  //    params.addRequiredCoupledVar("lambda", "Lagrange multiplier");
  //    params.addRequiredCoupledVar("conc", "concentration");
  params.addCoupledVar("cb", "cb");
  params.addCoupledVar("phi", "phi");
  params.addRequiredParam<MaterialPropertyName>("z", "The charge");
  params.addRequiredParam<Real>("kT", "The name of the thermal energy");
  return params;
}

ChemPotValueHC::ChemPotValueHC(const InputParameters & parameters)
  : ADKernel(parameters),
    //        _eta_name(_var.name()),
    // _mutarget(getADMaterialProperty<Real>("mutarget")),
    _mutarget(getParam<Real>("mutarget")),
    _mu0(getADMaterialProperty<Real>("mu0")),
    _cb(adCoupledValue("cb")),
    _phi(adCoupledValue("phi")),
    _z(getADMaterialProperty<Real>("z")),
    _kT(getParam<Real>("kT"))
//         _dh(getMaterialPropertyDerivative<Real>("h_name", _eta_name)),
//         _d2h(getMaterialPropertyDerivative<Real>("h_name", _eta_name, _eta_name)),
//         _d2ha(isCoupled("args") ? coupledComponents("args") :
//         coupledComponents("coupled_variables")), _d2ha_map(isCoupled("args") ?
//         getParameterJvarMap("args")
//                                                                 :
//                                                                 getParameterJvarMap("coupled_variables")),
//         _lambda(coupledValue("lambda")),
//         _lambda_var(coupled("lambda"))
{
  //     for (std::size_t i = 0; i < _d2ha.size(); ++i)
  //     {
  //         if (isCoupled("args"))
  //             _d2ha[i] = &getMaterialPropertyDerivative<Real>("h_name", _eta_name,
  //             coupledName("args", i));
  //         else
  //             _d2ha[i] = &getMaterialPropertyDerivative<Real>(
  //                     "h_name", _eta_name, coupledName("coupled_variables", i));
  //     }
}

ADReal
ChemPotValueHC::computeQpResidual()
{
  ADReal _mu =
      _mu0[_qp] + _kT * std::log(_u[_qp] / (1.0 - _u[_qp])) + _z[_qp] * 1.602e-19 * _phi[_qp];
  return (_mu - _mutarget) * _test[_i][_qp];
}

// Real
// ChemPotValueHC::computeQpJacobian()
// {
//     return _lambda[_qp] * _d2h[_qp] * _phi[_j][_qp] * _test[_i][_qp];
// }

// Real
// ChemPotValueHC::computeQpOffDiagJacobian(unsigned int jvar)
// {
//     if (jvar == _lambda_var)
//         return _phi[_j][_qp] * _dh[_qp] * _test[_i][_qp];

//     auto k = mapJvarToCvar(jvar, _d2ha_map);
//     if (k >= 0)
//         return _lambda[_qp] * (*_d2ha[k])[_qp] * _phi[_j][_qp] * _test[_i][_qp];

//     return 0.0;
// }
