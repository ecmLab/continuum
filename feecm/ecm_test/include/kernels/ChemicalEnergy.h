#ifndef CHEMICALENERGY_H
#define CHEMICALENERGY_H
#include "ADKernel.h"
class ChemicalEnergy : public ADKernel
{
public:
	static InputParameters validParams();
	ChemicalEnergy(const InputParameters & parameters);
protected:
	virtual ADReal computeQpResidual() override;
	const ADMaterialProperty<Real> & _w;
	const Real & _scale;
};
#endif //CHEMICALENERGY_H
