#pragma once
#include "ADKernel.h"
class FreeEnergyDouble : public ADKernel
{
public:
    static InputParameters validParams();
    FreeEnergyDouble(const InputParameters & parameters);
protected:
    virtual ADReal computeQpResidual() override;
    const ADMaterialProperty<Real> & _A;
    const Real & _scale;
};
