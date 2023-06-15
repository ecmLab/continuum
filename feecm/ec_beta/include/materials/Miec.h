
#pragma once

#include "ADMaterial.h"
/**
 * Material-derived objects override the computeQpProperties()
 * function.  They must declare and compute material properties for
 * use by other objects in the calculation such as Kernels and
 * BoundaryConditions.
 */
class Miec : public ADMaterial
{
public:
  static InputParameters validParams();

  Miec(const InputParameters & parameters);

protected:
  /// Necessary override. This is where the values of the properties are computed.
  virtual void computeQpProperties() override;

  /// Values from the input file
  const Real & _inIonicConductivity;
  const Real & _inElectronicConductivity;
  const Real & _inGbConductivity;

  /// The ionic conductivity of Li+ in SE
  ADMaterialProperty<Real> & _ionic_conductivity;
  ADMaterialProperty<Real> & _electronic_conductivity;
  ADMaterialProperty<Real> & _gb_conductivity;
  
//  usingMaterialMembers;
//  using ADMaterial<compute_stage>::_communicator;
};
