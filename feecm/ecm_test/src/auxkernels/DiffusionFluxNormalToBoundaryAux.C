

#include "DiffusionFluxNormalToBoundaryAux.h"
#include "Assembly.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("ecmApp", DiffusionFluxNormalToBoundaryAux);

// defineLegacyParams(DiffusionFluxNormalToBoundaryAux);

InputParameters
DiffusionFluxNormalToBoundaryAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Compute components of flux vector for diffusion problems "
                             "$(\\vv{J} = -D \\nabla C)$ normal to a boundary.");
  params.addRequiredCoupledVar("diffusion_variable", "The name of the variable");
  params.addRequiredParam<MaterialPropertyName>(
      "diffusivity",
      "The name of the diffusivity material property that will be used in the flux computation.");
  return params;
}

DiffusionFluxNormalToBoundaryAux::DiffusionFluxNormalToBoundaryAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _normals(_assembly.normals()),
    _grad_u(coupledGradient("diffusion_variable")),
    _diffusion_coef(getADMaterialProperty<Real>("diffusivity"))
{
}

Real
DiffusionFluxNormalToBoundaryAux::computeValue()
{
//    auto dim =
    RealVectorValue flux;
    for (int i =0; i < LIBMESH_DIM; i++)
    {
        flux(i) = MetaPhysicL::raw_value(-_diffusion_coef[_qp]) * _grad_u[_qp](i);
    }
  return flux * _normals[_qp];
}
