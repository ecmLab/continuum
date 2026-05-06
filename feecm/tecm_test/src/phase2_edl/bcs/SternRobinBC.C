#include "SternRobinBC.h"
#include "libmesh/vector_value.h"

registerADMooseObject("tecm_testApp", SternRobinBC);

InputParameters SternRobinBC::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  // Implement by substituting epsilon * dndu = C_S * (phi_M - u),
  // with C_S = epsilon0 * stern_relative_permittivity / stern_thickness
  params.addRequiredParam<Real>("epsilon0", "Vacuum permittivity epsilon_0 [F/m].");
  params.addRequiredParam<Real>("stern_relative_permittivity", "Relative permittivity of the Stern layer eps_S [-].");
  params.addRequiredParam<Real>("stern_thickness", "Stern layer thickness x_S [m].");
  params.addRequiredParam<Real>("metal_potential", "Electrode (metal) potential phi_M [V].");
  return params;
}

SternRobinBC::SternRobinBC(const InputParameters & parameters)
  : ADIntegratedBC(parameters),
    _C_S(getParam<Real>("epsilon0") *
         getParam<Real>("stern_relative_permittivity") /
         getParam<Real>("stern_thickness")),
    _phi_M(getParam<Real>("metal_potential"))
{
}

ADReal SternRobinBC::computeQpResidual()
{
  return -_C_S * (_phi_M - _u[_qp]) * _test[_i][_qp];
}
