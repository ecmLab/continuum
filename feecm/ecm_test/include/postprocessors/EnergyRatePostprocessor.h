
#pragma once

#include "GeneralPostprocessor.h"

 class EnergyRatePostprocessor : public GeneralPostprocessor
 {
 public:
   EnergyRatePostprocessor(const InputParameters & parameters);
  static InputParameters validParams();

   virtual void initialize();
   virtual void execute();
   virtual Real getValue() const;
 protected:
   const PostprocessorValue & _postprocessor, & _postprocessor_old, & _dt, & _dt_old;
 };

