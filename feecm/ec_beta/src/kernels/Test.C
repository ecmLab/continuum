#include "Test.h"

registerADMooseObject("ecBetaApp", Test);
InputParameters
Test::validParams()
{
    InputParameters params=ADKernel::validParams();
    params.addClassDescription("Free Energy Implimentation" );
    params.addParam<Real>("A",1,"a randomVal");
    return params;
}

Test::Test(const InputParameters & parameters) : ADKernel(parameters),
_A(getParam<Real>("A"))
{
}

ADReal
Test::computeQpResidual()
{
    return _A*_grad_test[_i][_qp]*_grad_u[_qp];
}
