
#include "NeumannEn.h"

registerADMooseObject("electro_chemo_mechApp", NeumannEn);

InputParameters
NeumannEn::validParams()
{
   
  InputParameters params = ADIntegratedBC::validParams();
  params.addClassDescription("The Neumann boundary condition for electron at Cathode side.");
  params.addRequiredParam<MaterialPropertyName>("LiPotElectrode", "The Li potental in anode or cathode, in V.");

// Add a coupled parameter: potLi
  params.addRequiredCoupledVar("potLi", "The variable representing the potential of Li+.");

// Add a parameter with a default value; this value can be overridden in the input file.
  params.addParam<Real>(
        "F_RT",
        38.68,
        "The constant of F/RT, when T = 300K.");
 
  return params;
}

NeumannEn::NeumannEn(const InputParameters & parameters)
  : ADIntegratedBC(parameters),

  // Couple to the potential of Li+
   _potLi(adCoupledValue("potLi")),
   _potLi_grad(adCoupledGradient("potLi")),

  // Get the parameters from the Material object
   _inlet_current(getADMaterialProperty<Real>("inlet_current")),
   _exchange_current(getADMaterialProperty<Real>("exchange_current")),
   _reaction_rate(getADMaterialProperty<Real>("reaction_rate")),

  // Get the parameters from the input file
   _LiPotElectrode(parameters.get<MaterialPropertyName>("LiPotElectrode")),
   _LiPotEle(getADMaterialProperty<Real>(_LiPotElectrode)),
   _F_RT(getParam<Real>("F_RT"))
{
}

ADReal
NeumannEn::computeQpResidual()
{
  ADReal k1 = std::exp(_reaction_rate[_qp] * _F_RT * (_potLi[_qp] + _u[_qp] - _LiPotEle[_qp]));
  ADReal k2 = std::exp(- (1 - _reaction_rate[_qp]) * _F_RT * (_potLi[_qp] + _u[_qp] - _LiPotEle[_qp]));

  return _test[_i][_qp] * (_inlet_current[_qp] + _exchange_current[_qp] * (k1 - k2));
}
