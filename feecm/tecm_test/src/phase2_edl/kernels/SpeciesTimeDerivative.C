#include "SpeciesTimeDerivative.h"

registerADMooseObject("tecm_testApp", SpeciesTimeDerivative);

InputParameters SpeciesTimeDerivative::validParams()
{
  InputParameters params = ADTimeDerivative::validParams();
  params.addClassDescription("Transient accumulation term for species concentration in the EDL model.");
  return params;
}

SpeciesTimeDerivative::SpeciesTimeDerivative(const InputParameters & parameters)
  : ADTimeDerivative(parameters)
{
}

ADReal SpeciesTimeDerivative::precomputeQpResidual()
{
  return ADTimeDerivative::precomputeQpResidual();
}

