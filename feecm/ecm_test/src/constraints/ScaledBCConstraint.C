
#include "ScaledBCConstraint.h"

registerADMooseObject("ecmApp", ScaledBCConstraint);

InputParameters
ScaledBCConstraint::validParams()
{
	InputParameters params = ADMortarConstraint::validParams();
    params.addClassDescription("Scales an existing contact problem variable"
                               "to apply on another variable");
    params.addRequiredParam<Real>("scale", "Scaling of flux");
    params.addParam<bool>("primary", true, "Apply bc on primary surface");
    params.set<bool>("compute_lm_residuals") = false;

    return params;
}

ScaledBCConstraint::ScaledBCConstraint(const InputParameters & parameters)
  : ADMortarConstraint(parameters), _scale(getParam<Real>("scale")),
        _primary(getParam<bool>("primary"))
{
}

ADReal
ScaledBCConstraint::computeQpResidual(Moose::MortarType mortar_type)
{
  switch (mortar_type)
  {
    case Moose::MortarType::Primary:
    {
        if (_primary)
            return -_lambda[_qp] * _test_primary[_i][_qp] * _scale;
    }
    case Moose::MortarType::Secondary:
    {
        if (!_primary)
            return -_lambda[_qp] * _test_secondary[_i][_qp] * _scale;
    }
    default:
      return 0;
  }
}
