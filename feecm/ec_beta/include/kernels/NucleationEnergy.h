#ifndef NUCLEATIONENERGY_H
#define NUCLEATIONENERGY_H

#include "ADKernel.h"

class NucleationEnergy : public ADKernel 
{
public:
	static InputParameters validParams();
	NucleationEnergy(const InputParameters & parameters);
protected:
	virtual ADReal computeQpResidual() override;
	const ADMaterialProperty<Real> & _delta_G;
	const ADMaterialProperty<Real> & _gamma;
	const Real & _scale;
};
#endif
