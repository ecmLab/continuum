#include "ADMaterial.h"
//#include "RealTensorValue.h"
#include "RotationTensor.h"

#pragma once
class RotatedConductivityTensor : public Material 
{
    public:
        static InputParameters validParams();
        RotatedConductivityTensor(const InputParameters & parameters);
    protected:
        virtual void computeQpProperties() override;
    private:
        RealTensorValue _conductivity_tensor;
        const Real & _sig_xx;
        const Real & _sig_xy;
        const Real & _sig_xz;
        const Real & _sig_yx;
        const Real & _sig_yy;
        const Real & _sig_yz;
        const Real & _sig_zx;
        const Real & _sig_zy;
        const Real & _sig_zz;
        const Real & _euler_angle_1;
        const Real & _euler_angle_2;
        const Real & _euler_angle_3;
        ADMaterialProperty<RealTensorValue> & _rotated_conductivity;
        ADMaterialProperty<Real> & _rotated_sig_xx;
        ADMaterialProperty<Real> & _rotated_sig_xy;
        ADMaterialProperty<Real> & _rotated_sig_xz;
        ADMaterialProperty<Real> & _rotated_sig_yx;
        ADMaterialProperty<Real> & _rotated_sig_yy;
        ADMaterialProperty<Real> & _rotated_sig_yz;
        ADMaterialProperty<Real> & _rotated_sig_zx;
        ADMaterialProperty<Real> & _rotated_sig_zy;
        ADMaterialProperty<Real> & _rotated_sig_zz;

};