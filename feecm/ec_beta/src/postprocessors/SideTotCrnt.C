#include "SideTotCrnt.h"

registerMooseObject("ecBetaApp", SideTotCrnt);

InputParameters
SideTotCrnt::validParams()
{
  InputParameters params = SideIntegralVariablePostprocessor::validParams();
  params.addRequiredParam<MaterialPropertyName>(
      "conductivity",
      "The name of the diffusivity material property that will be used in the flux computation.");
  params.addClassDescription("Computes the integral of the flux over the specified boundary");
  return params;
}

SideTotCrnt::SideTotCrnt(const InputParameters & parameters)
  : SideIntegralVariablePostprocessor(parameters),
    _conductivity(parameters.get<MaterialPropertyName>("conductivity")),
    _conductivity_coef(getADMaterialProperty<Real>(_conductivity))
{
}

Real
// Unit of total current in 2D: nA/um, in 3D: nA
SideTotCrnt::computeQpIntegral()
{

  Real k = 100*MetaPhysicL::raw_value(_conductivity_coef[_qp]) * _grad_u[_qp] * _normals[_qp];

  if(_conductivity == "ionic_conductivity") {
    k = -k;
  }

  return k;

}
