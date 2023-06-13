
#include "DirichletLi.h"

registerADMooseObject("ecBetaApp", DirichletLi);

InputParameters
DirichletLi::validParams()
{
  InputParameters params = ADNodalBC::validParams();

   // Add a coupled parameter: potEn
   params.addRequiredCoupledVar("potEn", "The potential of Electron");

  return params;
}

DirichletLi::DirichletLi(const InputParameters & parameters)
  : ADNodalBC(parameters),

    // Couple to the potential of Electron
     _potEn(adCoupledValue("potEn"))
{
}

ADReal
DirichletLi::computeQpResidual()
{
  // For dirichlet BCS, u=BC at the boundary, so the residual includes _u and the desired BC value:
//  return _u[_qp] - _potEn[_qp];
  return _u - _potEn[0];
}
