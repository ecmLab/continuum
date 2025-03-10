#ifndef COUPLINGFREEENERGY_H
#define COUPLINGFREEENERGY_H
#include "ADKernel.h"
class CouplingFreeEnergy : public ADKernel
{
public:
	static InputParameters validParams();
	CouplingFreeEnergy(const InputParameters & parameters);
protected:
	virtual ADReal computeQpResidual() override;
	const ADMaterialProperty<Real> & _c;
	const Real & _scale;
	const VariableValue & _coupledVar;
};
#endif //COUPLINGFREEENERGY_H
