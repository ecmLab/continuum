#include "FreeEnergyDouble.h"

registerADMooseObject("ecBetaApp", FreeEnergyDouble);
InputParameters
FreeEnergyDouble::validParams()
{
    InputParameters params=ADKernel::validParams();
    params.addClassDescription("Free Energy Implimentation" );
    params.addParam<MaterialPropertyName>("A",1,"Depth of The Double Well");
    params.addParam<Real>("scale",1,"Scaling Factor");
    return params;
}

FreeEnergyDouble::FreeEnergyDouble(const InputParameters & parameters) : ADKernel(parameters),
_A(getADMaterialProperty<Real>("A")),
_scale(getParam<Real>("scale"))
{
}

ADReal
FreeEnergyDouble::computeQpResidual()
{
    return _scale*_test[_i][_qp]*(2*_A[_qp]*std::pow(1 - _u[_qp],2)*_u[_qp] - 2*_A[_qp]*(1 - _u[_qp])*std::pow(_u[_qp],2));

}
