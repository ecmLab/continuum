#include "RotatedConductivityTensor.h"

registerMooseObject("ecBetaApp", RotatedConductivityTensor);
InputParameters
RotatedConductivityTensor::validParams()
{
    InputParameters params = Material::validParams();
    params.addClassDescription("Computes a Rotated conductivity Tensor");
    params.addRequiredParam<Real>("sig_xx", "the conductivity in xx direction");
    params.addRequiredParam<Real>("sig_xy", "the conductivity in xy direction");
    params.addRequiredParam<Real>("sig_xz", "the conductivity in xz direction");
    params.addRequiredParam<Real>("sig_yx", "the conductivity in yx direction");
    params.addRequiredParam<Real>("sig_yy", "the conductivity in yy direction");
    params.addRequiredParam<Real>("sig_yz", "the conductivity in yz direction");
    params.addRequiredParam<Real>("sig_zx", "the conductivity in zx direction");
    params.addRequiredParam<Real>("sig_zy", "the conductivity in zy direction");
    params.addRequiredParam<Real>("sig_zz", "the conductivity in zz direction");
    params.addRequiredParam<Real>("euler_angle_1", "First Euler angle (in radians)");
    params.addRequiredParam<Real>("euler_angle_2", "Second Euler angle (in radians)");
    params.addRequiredParam<Real>("euler_angle_3", "Third Euler angle (in radians)");
    return params;
}
RotatedConductivityTensor::RotatedConductivityTensor(const InputParameters & parameters)
: Material(parameters),
_sig_xx(getParam<Real>("sig_xx")),
_sig_xy(getParam<Real>("sig_xy")),
_sig_xz(getParam<Real>("sig_xz")),
_sig_yx(getParam<Real>("sig_yx")),
_sig_yy(getParam<Real>("sig_yy")),
_sig_yz(getParam<Real>("sig_yz")),
_sig_zx(getParam<Real>("sig_zx")),
_sig_zy(getParam<Real>("sig_zy")),
_sig_zz(getParam<Real>("sig_zz")),
_euler_angle_1(getParam<Real>("euler_angle_1")),
_euler_angle_2(getParam<Real>("euler_angle_2")),
_euler_angle_3(getParam<Real>("euler_angle_3")),
_rotated_conductivity(declareADProperty<RealTensorValue>("rotated_conductivity")),
_rotated_sig_xx(declareADProperty<Real>("rotated_sig_xx")),
_rotated_sig_xy(declareADProperty<Real>("rotated_sig_xy")),
_rotated_sig_xz(declareADProperty<Real>("rotated_sig_xz")),
_rotated_sig_yx(declareADProperty<Real>("rotated_sig_yx")),
_rotated_sig_yy(declareADProperty<Real>("rotated_sig_yy")),
_rotated_sig_yz(declareADProperty<Real>("rotated_sig_yz")),
_rotated_sig_zx(declareADProperty<Real>("rotated_sig_zx")),
_rotated_sig_zy(declareADProperty<Real>("rotated_sig_zy")),
_rotated_sig_zz(declareADProperty<Real>("rotated_sig_zz"))
{
    _conductivity_tensor = RealTensorValue(_sig_xx, _sig_xy, _sig_xz,
                                           _sig_yx, _sig_yy, _sig_yz, 
                                           _sig_zx, _sig_zy, _sig_zz);
}
void 
RotatedConductivityTensor::computeQpProperties()
{
    RealVectorValue euler_angles(_euler_angle_1, _euler_angle_2, _euler_angle_3);
    RotationTensor rotation(euler_angles);
    _rotated_conductivity[_qp] = rotation.transpose() * _conductivity_tensor * rotation;
    _rotated_sig_xx[_qp] = _rotated_conductivity[_qp](0,0);
    _rotated_sig_xy[_qp] = _rotated_conductivity[_qp](0,1);
    _rotated_sig_xz[_qp] = _rotated_conductivity[_qp](0,2);
    _rotated_sig_yx[_qp] = _rotated_conductivity[_qp](1,0);
    _rotated_sig_yy[_qp] = _rotated_conductivity[_qp](1,1);
    _rotated_sig_yz[_qp] = _rotated_conductivity[_qp](1,2);
    _rotated_sig_zx[_qp] = _rotated_conductivity[_qp](2,0);
    _rotated_sig_zy[_qp] = _rotated_conductivity[_qp](2,1);
    _rotated_sig_zz[_qp] = _rotated_conductivity[_qp](2,2);
}