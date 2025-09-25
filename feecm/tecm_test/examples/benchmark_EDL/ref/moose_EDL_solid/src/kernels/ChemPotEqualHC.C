// Created by: Zeeshan Ahmad

#include "ChemPotEqualHC.h"

registerMooseObject("edl_solidApp", ChemPotEqualHC);

InputParameters
ChemPotEqualHC::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription(
      "Kernel to enforce pointwise equilibrium of (electro)chemical potentials mua + mub = 0.");
  params.addRequiredParam<MaterialPropertyName>("mu0a", "The name of the first chemical potential");
  params.addRequiredParam<MaterialPropertyName>("mu0b",
                                                "The name of the second chemical potential");
  //    params.addRequiredCoupledVar("lambda", "Lagrange multiplier");
  //    params.addRequiredCoupledVar("conc", "concentration");
  params.addCoupledVar("cb", "cb");
  params.addCoupledVar("phi", "phi");
  params.addRequiredParam<MaterialPropertyName>("za", "The name of the first charge number");
  params.addRequiredParam<MaterialPropertyName>("zb", "The name of the second charge number");
  params.addRequiredParam<Real>("kT", "The name of the thermal energy");
  return params;
}

ChemPotEqualHC::ChemPotEqualHC(const InputParameters & parameters)
  : ADKernel(parameters),
    //        _eta_name(_var.name()),
    _mu0a(getADMaterialProperty<Real>("mu0a")),
    _mu0b(getADMaterialProperty<Real>("mu0b")),
    _cb(adCoupledValue("cb")),
    _phi(adCoupledValue("phi")),
    _za(getADMaterialProperty<Real>("za")),
    _zb(getADMaterialProperty<Real>("zb")),
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
ChemPotEqualHC::computeQpResidual()
{
  ADReal _mua =
      _mu0a[_qp] + _kT * std::log(_u[_qp] / (1.0 - _u[_qp])) + _za[_qp] * 1.602e-19 * _phi[_qp];
  ADReal _mub =
      _mu0b[_qp] + _kT * std::log(_cb[_qp] / (1.0 - _cb[_qp])) + _zb[_qp] * 1.602e-19 * _phi[_qp];
  return (_mua + _mub) * _test[_i][_qp];
}

// Real
// ChemPotEqualHC::computeQpJacobian()
// {
//     return _lambda[_qp] * _d2h[_qp] * _phi[_j][_qp] * _test[_i][_qp];
// }

// Real
// ChemPotEqualHC::computeQpOffDiagJacobian(unsigned int jvar)
// {
//     if (jvar == _lambda_var)
//         return _phi[_j][_qp] * _dh[_qp] * _test[_i][_qp];

//     auto k = mapJvarToCvar(jvar, _d2ha_map);
//     if (k >= 0)
//         return _lambda[_qp] * (*_d2ha[k])[_qp] * _phi[_j][_qp] * _test[_i][_qp];

//     return 0.0;
// }
