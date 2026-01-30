#pragma once

#include "ADTimeDerivative.h"

/**
 * Simple AD time-derivative kernel for species transport in the Phase 2 EDL model:
 * residual contribution = (∂c/∂t, w).
 */
class SpeciesTimeDerivative : public ADTimeDerivative
{
public:
  static InputParameters validParams();
  SpeciesTimeDerivative(const InputParameters & parameters);

protected:
  ADReal precomputeQpResidual() override;
};

