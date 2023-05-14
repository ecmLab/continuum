
#include "DirichletEn.h"

registerADMooseObject("electro_chemo_mechApp", DirichletEn);

InputParameters
DirichletEn::validParams()
{
  InputParameters params = ADNodalBC::validParams();

   // Add a coupled parameter: potLi
  params.addRequiredCoupledVar("potLi", "The potential of Li-ion");
  return params;
}

DirichletEn::DirichletEn(const InputParameters & parameters)
  : ADNodalBC(parameters),

    // Couple to the potential of Li+
     _potLi(adCoupledValue("potLi"))
{
}

ADReal
DirichletEn::computeQpResidual()
{
  // For dirichlet BCS, u=BC at the boundary, so the residual includes _u and the desired BC value:
//  return _u[_qp] - _potLi[_qp];
  return _u - _potLi[0];
}
