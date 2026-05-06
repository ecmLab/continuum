#include "ADNernstPlanckConvection.h"

#include "MooseUtils.h"

registerADMooseObject("tecm_testApp", ADNernstPlanckConvection);

InputParameters
ADNernstPlanckConvection::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Electromigration term for binary species transport in the EDL model.");

  params.addRequiredParam<MaterialPropertyName>("diffusivity", "Diffusivity material property D_i [m²/s].");
  params.addCoupledVar("potential", "Electric potential variable φ.");
  params.addRequiredParam<Real>("valence", "Integer valence z_i.");
  params.addRequiredParam<Real>("faraday_constant", "Faraday constant F in C/mol.");
  params.addRequiredParam<Real>("gas_constant", "Universal gas constant R.");
  params.addCoupledVar("temperature", "Temperature variable (K).");
  return params;
}

ADNernstPlanckConvection::ADNernstPlanckConvection(const InputParameters & parameters)
  : ADKernel(parameters),
    _diffusivity(getADMaterialProperty<Real>("diffusivity")),
    _grad_potential(adCoupledGradient("potential")),
    _temperature(adCoupledValue("temperature")),
    _valence(getParam<Real>("valence")),
    _faraday_constant(getParam<Real>("faraday_constant")),
    _gas_constant(getParam<Real>("gas_constant"))
{
  if (MooseUtils::absoluteFuzzyEqual(_gas_constant, 0.0))
    paramError("gas_constant", "Gas constant must be non-zero.");
}

ADReal
ADNernstPlanckConvection::computeQpResidual()
{
  const ADReal mobility = _valence * _faraday_constant * _diffusivity[_qp] /
                           (_gas_constant * _temperature[_qp]);

  const ADRealVectorValue flux = -mobility * _u[_qp] * _grad_potential[_qp];
  return flux * _grad_test[_i][_qp];
}
