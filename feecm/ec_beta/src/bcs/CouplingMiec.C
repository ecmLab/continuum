
#include "CouplingMiec.h"

registerADMooseObject("ecBetaApp", CouplingMiec);

InputParameters
CouplingMiec::validParams()
{
 InputParameters params = ADIntegratedBC::validParams();
 params.addClassDescription("Electron and Li-ions at interface and surface of the Miec are coupled through the B-V relation.");

// Add a parameter with a default value; this value can be overridden in the input file.
 params.addRequiredParam<Real>("LiCrtRef", "The reference current at the interface/surface, in unit mA/cm^2.");

// Add a coupled parameter: potLi
 params.addRequiredCoupledVar("potLi", "The potential of Li+ in SE");

 return params;
}

CouplingMiec::CouplingMiec(const InputParameters & parameters)
  : ADIntegratedBC(parameters),

  // Get the parameters from the input file
   _liCrtRef(getParam<Real>("LiCrtRef")),

  // Get the parameters from the material property
   _ionic_conductivity(getADMaterialProperty<Real>("ionic_conductivity")),

  // Couple to the Li potential
   _potLi(adCoupledValue("potLi")),
  // Get the gradient of the Li potential
   _potLi_gradient(adCoupledGradient("potLi"))
{
}

ADReal
CouplingMiec::computeQpResidual()
{
  return -_test[_i][_qp] * (_liCrtRef + 10*_ionic_conductivity[_qp] * _potLi_gradient[_qp] * _normals[_qp]);

}
