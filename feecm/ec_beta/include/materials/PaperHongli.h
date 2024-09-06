#pragma once

#include "Material.h"

/**
 * Material objects inherit from Material and override computeQpProperties.
 *
 * Their job is to declare properties for use by other objects in the
 * calculation such as Kernels and BoundaryConditions.
 */
class PaperHongli : public Material
{
public:
  static InputParameters validParams();

  PaperHongli(const InputParameters & parameters);

protected:
  /// Necessary override. This is where the values of the properties are computed.
  virtual void computeQpProperties() override;

  const ADVariableValue & _cLi;
  const Real & _c1; 
  const Real & _c2;

  /// The diffusivity um^2/s
  ADMaterialProperty<Real> & _diffusivity;
  /// The conductivity mS/cm
  ADMaterialProperty<Real> & _ionic_conductivity;

};
