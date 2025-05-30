#pragma once

#include "ADMaterial.h"
#include "RankTwoTensor.h"
//#include "ADRankTwoTensor.h"

class AnisotropyMaterialConst  : public ADMaterial
{
public:
    static InputParameters validParams();
    AnisotropyMaterialConst(const InputParameters & parameters);

protected:
    virtual void computeQpProperties() override;

private:
  
  const Real & _sig11;
  const Real & _sig12;
  const Real & _sig13;
  const Real & _sig21;
  const Real & _sig22;
  const Real & _sig23;
  const Real & _sig31;
  const Real & _sig32;
  const Real & _sig33;
  ADMaterialProperty<RealTensorValue> & _tensor_property;


};
//InputParameters validParams<AnisotropyMaterialConst<T>>();
