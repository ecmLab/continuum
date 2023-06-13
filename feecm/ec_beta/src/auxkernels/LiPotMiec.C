
#include "LiPotMiec.h"

registerMooseObject("ecBetaApp", LiPotMiec);

InputParameters
LiPotMiec::validParams()
{
  InputParameters params = AuxKernel::validParams();

  // Add a "coupling paramater" to get a variable from the input file.
  params.addRequiredCoupledVar("potLi", "The potential field of Li-ion.");

  // Add a "coupling paramater" to get a variable from the input file.
  params.addRequiredCoupledVar("potEn", "The potential field of Electron.");

  return params;
}

LiPotMiec::LiPotMiec(const InputParameters & parameters)
  : AuxKernel(parameters),

    // Couple to the potential of Li+
   _potLi(coupledValue("potLi")),

    // Couple to the potential of En
   _potEn(coupledValue("potEn"))

{
}

Real
LiPotMiec::computeValue()
{
  Real kk = _potLi[_qp] + _potEn[_qp];
//  if (kk<0) {
//    kk = 0.0;
//  } 
  return kk;
}
