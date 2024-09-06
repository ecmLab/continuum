
#pragma once

// MOOSE includes
#include "SideIntegralVariablePostprocessor.h"

// Forward Declarations
/**
 * This postprocessor computes a side integral of the mass flux.
 */
class SideTotCrnt : public SideIntegralVariablePostprocessor
{
public:
  static InputParameters validParams();

  SideTotCrnt(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  MaterialPropertyName _conductivity;
  const ADMaterialProperty<Real> & _conductivity_coef;
};

