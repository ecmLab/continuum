
#pragma once

#include "ADMaterial.h"
/**
 * Material-derived objects override the computeQpProperties()
 * function.  They must declare and compute material properties for
 * use by other objects in the calculation such as Kernels and
 * BoundaryConditions.
 */
class Ionics : public ADMaterial
{
public:
  static InputParameters validParams();

  Ionics(const InputParameters & parameters);

protected:
  /// Necessary override. This is where the values of the properties are computed.
  virtual void computeQpProperties() override;

  /// Values from the input file
  const Real & _inIonicConductivity;
  const Real & _inGbConductivity;
//  const Real & _inExchangeCurrent;
//  const Real & _inReactionRate;
//  const Real & _exchange_current;

  /// The ionic conductivity of Li+ in SE
  ADMaterialProperty<Real> & _ionic_conductivity;
  ADMaterialProperty<Real> & _gb_conductivity;
  
  /// The exchange current density at Li metal/SE interface
//  ADMaterialProperty<Real> & _exchange_current;
  
  /// The reaction rate of the Li/Li+ reaction
//  ADMaterialProperty<Real> & _reaction_rate;


//  usingMaterialMembers;
//  using ADMaterial<compute_stage>::_communicator;
};
