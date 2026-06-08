#pragma once
#include "Material.h"
#include "DerivativeMaterialInterface.h"
class ChemicalEnergyDensity_1 : public DerivativeMaterialInterface<Material>
{
public:
static InputParameters validParams();
ChemicalEnergyDensity_1(const InputParameters & parameters);
protected:
const MaterialPropertyName _energy_name;
const MooseVariable *_c_var;
const ADVariableValue & _c;
const ADVariableValue & _c_dot;
 /// The chemical energy density rate
  ADMaterialProperty<Real> & _psi_dot;

  /// Derivative of the chemical energy density rate w.r.t. the chemical concentration rate
  ADMaterialProperty<Real> & _d_psi_dot_d_c_dot;
  };
