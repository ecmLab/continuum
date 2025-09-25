#include "ADNernstPlanckConvection.h"

registerMooseObject("tecm_testApp", ADNernstPlanckConvection);

InputParameters
ADNernstPlanckConvection::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Nernst-Planck convection kernel: v · ∇c where v = -μ * z * ∇Ψ");
  
  params.addRequiredParam<MaterialPropertyName>("diffusivity", "The diffusivity material property");
  params.addRequiredParam<Real>("valence", "Ion valence (charge number)");
  params.addRequiredParam<Real>("faraday_constant", "Faraday constant [C/mol]");
  params.addRequiredParam<Real>("gas_constant", "Gas constant [J/(mol·K)]");
  params.addRequiredParam<Real>("temperature", "Temperature [K]");
  params.addRequiredCoupledVar("potential", "Electric potential variable");
  
  return params;
}

ADNernstPlanckConvection::ADNernstPlanckConvection(const InputParameters & parameters)
  : ADKernel(parameters),
    _diffusivity(getADMaterialProperty<Real>("diffusivity")),
    _valence(getParam<Real>("valence")),
    _faraday_constant(getParam<Real>("faraday_constant")),
    _gas_constant(getParam<Real>("gas_constant")),
    _temperature(getParam<Real>("temperature")),
    _potential_gradient(adCoupledGradient("potential"))
{
}

ADReal
ADNernstPlanckConvection::computeQpResidual()
{
  // Electrophoretic mobility: μ = D * z * F / (RT) [m²/(V·s)]
  ADReal mobility = _diffusivity[_qp] * _valence * _faraday_constant / (_gas_constant * _temperature);
  
  // Ion velocity: v = -μ * ∇Ψ [m/s]
  // Negative sign because positive ions move toward lower potential
  ADRealVectorValue ion_velocity = -mobility * _potential_gradient[_qp];
  
  // Convection term: v · ∇c
  // Weak form: ∫ (v · ∇u) * test dV
  // This gives the correct electromigration transport
  return ion_velocity * _grad_u[_qp] * _test[_i][_qp];
}