
#include "EnergyRatePostprocessor.h"

registerMooseObject("ecBetaApp", EnergyRatePostprocessor);

InputParameters EnergyRatePostprocessor::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addClassDescription("Calculates the change of a postprocessor divided by the time step.");
  params.addRequiredParam<PostprocessorName>("postprocessor", "The name of the postprocessor used for exit criterion");
  params.addRequiredParam<PostprocessorName>("dt", "The dt postprocessor");
  return params;
}

EnergyRatePostprocessor::EnergyRatePostprocessor(const InputParameters & parameters) :
    GeneralPostprocessor(parameters),
    _postprocessor(getPostprocessorValue("postprocessor")),
    _postprocessor_old(getPostprocessorValueOld("postprocessor")),
    _dt(getPostprocessorValue("dt")),
    _dt_old(getPostprocessorValueOld("dt"))
{
}

void
EnergyRatePostprocessor::initialize(){
}

void
EnergyRatePostprocessor::execute(){
}

Real
EnergyRatePostprocessor::getValue() const
{
  return fabs( ( fabs(_postprocessor)- fabs(_postprocessor_old) )*pow(fabs(_postprocessor),-1)) / _dt;
}
