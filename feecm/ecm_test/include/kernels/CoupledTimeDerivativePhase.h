#ifndef COUPLEDTIMEDERIVATIVEPHASE_H
#define COUPLEDTIMEDERIVATIVEPHASE_H
#include "ADKernel.h"

class CoupledTimeDerivativePhase : public ADKernel
{
public:
	static InputParameters validParams();
	CoupledTimeDerivativePhase(const InputParameters & parameters);
protected:
	virtual ADReal computeQpResidual() override;
	const Real & _n;
	const Real & _F;
	const Real & _scale;
	const ADMaterialProperty<Real> & _C;
	const ADVariableValue & _phi_var_dot; 
};
#endif
