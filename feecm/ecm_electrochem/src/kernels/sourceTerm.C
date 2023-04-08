#include "sourceTerm.h"

#include "Function.h"

registerADMooseObject("liExpulsionApp", sourceTerm);

InputParameters
sourceTerm::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Implements the weak form $(\\psi_i, -f)$.");

// Add a coupled parameter: Concentration of the ions
  params.addRequiredParam<MaterialPropertyName>("srcCoef", "The source coefficient.");
  params.addParam<FunctionName>("function", "0", "A function that describes the source");
  params.addParam<Real>("F_RT", 38.68, "The constant of F/RT,in unit 1/V, when T = 300K.");
    
  return params;
}

sourceTerm::sourceTerm(const InputParameters & parameters)
  : Kernel(parameters),
   _srcCoef(getMaterialProperty<Real>("srcCoef")),
   _function(getFunction("function")),
   _F_RT(getParam<Real>("F_RT"))
{
}

Real
sourceTerm::computeQpResidual()
{
  return _test[_i][_qp] * (_F_RT * _srcCoef[_qp] * _u[_qp] * _u[_qp] - _function.value(_t, _q_point[_qp]));
}

Real
sourceTerm::computeQpJacobian()
{
  return _test[_i][_qp] * _F_RT * _srcCoef[_qp] * 2*_u[_qp] * _phi[_j][_qp] ;
}
