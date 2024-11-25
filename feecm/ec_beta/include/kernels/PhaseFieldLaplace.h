#pragma once
#include "ADKernel.h"
class PhaseFieldLaplace : public ADKernel 
{
    public:
        static InputParameters validParams();
        PhaseFieldLaplace(const InputParameters & parameters);
    protected:
        virtual ADReal computeQpResidual() override;
        
        const ADMaterialProperty<Real> & _k0;
        
};
