#include "PoissonEquation.h"

registerADMooseObject("tecm_testApp", PoissonEquation);

InputParameters
PoissonEquation::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription(
      "Implements the space charge Poisson equation for TECM without electroneutrality: "
      "-∇ · (ε ∇Ψ) = -ρₑ, where ρₑ = F(z₊c₊ + z₋c₋) + ρfixed");
  
  params.addRequiredParam<MaterialPropertyName>("permittivity", 
                                                "Permittivity material property");
  params.addRequiredCoupledVar("charge_density", 
                               "Charge density variable (ρₑ)");
  
  return params;
}

PoissonEquation::PoissonEquation(const InputParameters & parameters)
  : ADKernel(parameters),
    _permittivity(getADMaterialProperty<Real>("permittivity")),
    _charge_density(adCoupledValue("charge_density"))
{
}

ADReal
PoissonEquation::computeQpResidual()
{
  // Implement: -∇ · (ε ∇Ψ) = -ρₑ
  // In weak form: ε ∇Ψ · ∇φ = ρₑ φ
  // This gives us: ε ∇u · ∇test - ρₑ test
  return _permittivity[_qp] * _grad_u[_qp] * _grad_test[_i][_qp] - 
         _charge_density[_qp] * _test[_i][_qp];
}