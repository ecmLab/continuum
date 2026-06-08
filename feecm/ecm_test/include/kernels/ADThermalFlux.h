#pragma once
#include "ADKernel.h"
class ADThermalFlux : public ADKernel 
{
    public:
        static InputParameters validParams();
        ADThermalFlux(const InputParameters & parameters);
    protected:
        virtual ADReal computeQpResidual() override;
        
        const MaterialProperty<Real> & _k;
        const Real & _scale;
        
};
