#include "NonDimensionalParameters.h"

registerMooseObject("tecm_testApp", NonDimensionalParameters);

InputParameters
NonDimensionalParameters::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material class providing characteristic scales and dimensionless "
                            "parameters for TECM non-dimensionalization framework");

  // Characteristic scale parameters (with typical Li-ion battery values as defaults)
  params.addParam<Real>("c0", 29400.0, "Reference concentration [mol/m³] (Li in Li metal)");
  params.addParam<Real>("T0", 298.0, "Reference temperature [K]");
  params.addParam<Real>("D0", 1.0e-14, "Reference diffusivity [m²/s] (Li in graphite)");
  params.addParam<Real>("j0", 1.0, "Exchange current density [A/m²] (Li metal/SE interface)");
  params.addParam<Real>("Omega0", 1.304e-5, "Li molar volume in Li metal [m³/mol]");
  params.addParam<Real>("F", 96485.0, "Faraday constant [C/mol]");
  params.addParam<Real>("R", 8.314, "Gas constant [J/mol/K]");

  // Optional material properties for local calculations
  params.addParam<MaterialPropertyName>("permittivity", "Local permittivity for electrostatic parameter");

  return params;
}

NonDimensionalParameters::NonDimensionalParameters(const InputParameters & parameters)
  : Material(parameters),
    // Input parameters
    _c0(getParam<Real>("c0")),
    _T0(getParam<Real>("T0")),
    _D0(getParam<Real>("D0")),
    _j0(getParam<Real>("j0")),
    _Omega0(getParam<Real>("Omega0")),
    _F(getParam<Real>("F")),
    _R(getParam<Real>("R")),
    // Characteristic scales
    _L0(declareProperty<Real>("L0")),
    _t0(declareProperty<Real>("t0")),
    _phi0(declareProperty<Real>("phi0")),
    _mu0(declareProperty<Real>("mu0")),
    _sigma0(declareProperty<Real>("sigma0")),
    // Dimensionless parameters
    _kappa(declareProperty<Real>("kappa")),
    _M_mech(declareProperty<Real>("M_mech")),
    // Optional properties
    _epsilon(isParamValid("permittivity") ? &getMaterialProperty<Real>("permittivity") : nullptr)
{
}

void
NonDimensionalParameters::computeQpProperties()
{
  // Characteristic length: L₀ = Fc₀D₀/j₀
  _L0[_qp] = _F * _c0 * _D0 / _j0;
  
  // Characteristic time: t₀ = L₀²/D₀ = F²c₀²D₀/j₀²
  _t0[_qp] = _L0[_qp] * _L0[_qp] / _D0;
  
  // Characteristic potential: φ₀ = RT₀/F (thermal voltage)
  _phi0[_qp] = _R * _T0 / _F;
  
  // Characteristic chemical potential: μ₀ = RT₀
  _mu0[_qp] = _R * _T0;
  
  // Characteristic stress: σ₀ = RT₀/Ω₀
  _sigma0[_qp] = _R * _T0 / _Omega0;
  
  // Mechanical coupling parameter: M = Ω₀c₀
  _M_mech[_qp] = _Omega0 * _c0;
  
  // Electrostatic parameter: κ = F²c₀L₀²/(εRT₀)
  if (_epsilon)
  {
    _kappa[_qp] = _F * _F * _c0 * _L0[_qp] * _L0[_qp] / ((*_epsilon)[_qp] * _R * _T0);
  }
  else
  {
    // Use vacuum permittivity as default: ε₀ = 8.854e-12 F/m
    const Real epsilon0 = 8.854e-12;
    _kappa[_qp] = _F * _F * _c0 * _L0[_qp] * _L0[_qp] / (epsilon0 * _R * _T0);
  }
}