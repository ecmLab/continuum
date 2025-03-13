#include "ChemicalEnergyDensity_1.h"


registerMooseObject("EelApp", ChemicalEnergyDensity_1);

InputParameters
ChemicalEnergyDensity_1::validParams()
{

  InputParameters params= DerivativeMaterialInterface<Material>::validParams();
  params.addRequiredCoupledVar("concentration","chemical_potential");
  params.addRequiredParam<MaterialPropertyName>("chemical_energy_density","Name of Chemical Energy Density");
  return params;
}
ChemicalEnergyDensity_1::ChemicalEnergyDensity_1(const InputParameters & parameters)
: DerivativeMaterialInterface<Material>(parameters),
_energy_name(getParam<MaterialPropertyName>("chemical_energy_density")),
_c_var(getVar("concentration", 0)),
_c(adCoupledValue("concentration")),
_c_dot(adCoupledDot("concentration")),
_psi_dot(declareADProperty<Real>("dot(" + _energy_name + ")")),
_d_psi_dot_d_c_dot(declarePropertyDerivative<Real, true>("dot(" + _energy_name + ")",
                                                         "dot(" + _c_var->name() + ")"))
{
}
