#ifndef ELECTRODEDRIVINGFORCE_H
#define ELECTRODEDRIVINGFORCE_H
#include "ADKernel.h"
class ElectrodeDrivingForce : public ADKernel
{
public:
	static InputParameters validParams();
	ElectrodeDrivingForce(const InputParameters & parameters);
protected:
	virtual ADReal computeQpResidual() override;
	const VariableValue & _pot;
	const VariableValue & _ref_pot;
	const VariableValue & _conc;
	const ADMaterialProperty<Real> & _h;
	const Real & _alpha;
	const Real & _beta;
	const Real & _n;
	const Real & _R;
	const Real & _T;
	const Real & _F;
	const Real & _scale;
	
};
#endif //ELECTRODEDRIVINGFORCE_H
