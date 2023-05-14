
#pragma once

// MOOSE includes
#include "SideIntegralVariablePostprocessor.h"

// Forward Declarations
class sideTotCrnt;

/**
 * This postprocessor computes a side integral of the mass flux.
 */
class sideTotCrnt : public SideIntegralVariablePostprocessor
{
public:
  static InputParameters validParams();

  sideTotCrnt(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  MaterialPropertyName _conductivity;
  const ADMaterialProperty<Real> & _conductivity_coef;
};

