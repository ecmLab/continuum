#include "EDLFluxes.h"

registerADMooseObject("tecm_testApp", EDLFluxes);

InputParameters EDLFluxes::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription("Provides fluxes for Poisson and Nernst-Planck equations in EDL.");
  params.addRequiredCoupledVar("phi", "Electric potential");
  params.addRequiredCoupledVar("c_plus", "Cation concentration");
  params.addRequiredCoupledVar("c_minus", "Anion concentration");

  params.addRequiredParam<MaterialPropertyName>("permittivity", "Permittivity property name ε");
  params.addRequiredParam<MaterialPropertyName>("D_plus", "Diffusivity of cation D_plus");
  params.addRequiredParam<MaterialPropertyName>("D_minus", "Diffusivity of anion D_minus");

  params.addParam<Real>("faraday_constant", 96485.33212, "Faraday constant F [C/mol]");
  params.addParam<Real>("gas_constant", 8.314462618, "Gas constant R [J/(mol·K)]");
  params.addParam<Real>("temperature", 298.15, "Absolute temperature T [K]");
  params.addParam<Real>("z_plus", 1.0, "Valence of cation");
  params.addParam<Real>("z_minus", -1.0, "Valence of anion");

  return params;
}

EDLFluxes::EDLFluxes(const InputParameters & parameters)
  : ADMaterial(parameters),
    _grad_phi(adCoupledGradient("phi")),
    _grad_c_plus(adCoupledGradient("c_plus")),
    _grad_c_minus(adCoupledGradient("c_minus")),
    _c_plus_val(adCoupledValue("c_plus")),
    _c_minus_val(adCoupledValue("c_minus")),
    _permittivity(getADMaterialProperty<Real>("permittivity")),
    _D_plus(getADMaterialProperty<Real>("D_plus")),
    _D_minus(getADMaterialProperty<Real>("D_minus")),
    _F(getParam<Real>("faraday_constant")),
    _R(getParam<Real>("gas_constant")),
    _T(getParam<Real>("temperature")),
    _z_plus(getParam<Real>("z_plus")),
    _z_minus(getParam<Real>("z_minus")),
    _elec_flux(declareADProperty<ADRealVectorValue>("elec_flux")),
    _j_diff_plus(declareADProperty<ADRealVectorValue>("j_diff_plus")),
    _j_diff_minus(declareADProperty<ADRealVectorValue>("j_diff_minus")),
    _j_em_plus(declareADProperty<ADRealVectorValue>("j_em_plus")),
    _j_em_minus(declareADProperty<ADRealVectorValue>("j_em_minus")),
    _charge_density(declareADProperty<ADReal>("charge_density"))
{
}

void EDLFluxes::computeQpProperties()
{
  // Electrical flux ε ∇φ
  _elec_flux[_qp] = _permittivity[_qp] * _grad_phi[_qp];

  // Diffusion fluxes
  _j_diff_plus[_qp] = -_D_plus[_qp] * _grad_c_plus[_qp];
  _j_diff_minus[_qp] = -_D_minus[_qp] * _grad_c_minus[_qp];

  // Electromigration fluxes: (z F D)/(R T) c ∇φ
  const ADReal coeff_p = (_z_plus * _F) * _D_plus[_qp] / (_R * _T);
  const ADReal coeff_m = (_z_minus * _F) * _D_minus[_qp] / (_R * _T);
  _j_em_plus[_qp] = -coeff_p * _c_plus_val[_qp] * _grad_phi[_qp];
  _j_em_minus[_qp] = -coeff_m * _c_minus_val[_qp] * _grad_phi[_qp];

  // Charge density
  _charge_density[_qp] = _F * (_z_plus * _c_plus_val[_qp] + _z_minus * _c_minus_val[_qp]);
}
