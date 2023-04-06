/*
 * PandolfiOrtizReversible.C
 *
 *  Created on: Jun 17, 2020
 *      Author: srinath
 */

#include "PandolfiOrtizReversible.h"

registerADMooseObject("electro_chemo_mechApp", PandolfiOrtizReversible);

InputParameters
PandolfiOrtizReversible::validParams()
{
	InputParameters params = CZMMaterialBase::validParams();
	params.addClassDescription("Ortiz and Pandolfi bi-linear cohesive law with damage");
	params.addRequiredParam<Real>("gap_at_max_normal_traction",
			"interface displacement jump when traction reaches max stress");
	params.addRequiredParam<Real>("gap_at_zero_normal_traction",
				"interface displacement jump when traction reaches zero");
	params.addRequiredParam<Real>("compressive_stiffness",
			"Compressive stiffness in normal opening mode to prevent crack faces overlap");
	params.addRequiredParam<Real>("undamaged_stiffness",
			"Tensile stiffness of undamaged interface");
	params.addRequiredParam<Real>("relative_stiffness",
			"Relative stiffness between shear and tension");
	params.addRequiredParam<Real>("viscosity",
			"Artificial viscosity for numerical stability");
	return params;
}

PandolfiOrtizReversible::PandolfiOrtizReversible(const InputParameters & parameters)
	: CZMMaterialBase(parameters),
	_d1(getParam<Real>("gap_at_max_normal_tracion")),
	_d2(getParam<Real>("gap_at_zero_normal_tracion")),
	_k1(getParam<Real>("compressive_stiffness")),
	_k0(getParam<Real>("undamaged_stiffness")),
	_theta(getParam<Real>("relative_stiffness")),
	_zeta(getParam<Real>("viscosity")),
	_sigma_max(_k0*_d2),
	_delta_eff(declareProperty<Real>("delta_eff")),
	_delta_eff_old(getMaterialPropertyOld<Real>("delta_eff")),
	_T_eff(declareProperty<Real>("T_eff")),
	_T_eff_old(getMaterialPropertyOld<Real>("T_eff")),
	_damage(declareProperty<Real>("damage")),
	_damage_old(getMaterialPropertyOld<Real>("damage")),
	_damage_strength(declareProperty<Real>("damage_stregnth")),
	_damage_strength_old(getMaterialPropertyOld<Real>("damage_strength"))
{
}

void
PandolfiOrtizReversible::initQpStatefulProperties()
{
	_delta_eff[_qp] = 0.0;
	_T_eff[_qp] = 0.0;
	_damage[_qp] = 0.0;
	_damage_strength[_qp] = _sigma_max;
}

RealVectorValue
PandolfiOrtizReversible::computeTraction()
{
	RealVectorValue traction_local;
	RealVectorValue displacement_jump_vel;
	displacement_jump_vel = (_displacement_jump[_qp] - _displacement_jump_old[_qp])/_dt;

	_damage_strength[_qp] = _sigma_max * (1.0 - _damage_old[_qp]) * _d2;
	_damage_strength[_qp] /= _d2 - _damage_old[_qp] * (_d2 - _d1);

	Real delta_eff;
	delta_eff = std::pow(_displacement_jump[_qp](0),2);
	delta_eff += std::pow(_theta,2) * std::pow(_displacement_jump[_qp](1),2);
	delta_eff += std::pow(_theta,2) * std::pow(_displacement_jump[_qp](2),2);
	_delta_eff[_qp] = std::sqrt(delta_eff);


	for (int i = 0; i < 3; i++)
	{
		if (i == 0)
		{
			traction_local(i) = _k0 * (1.0 -_damage_old[_qp]) * _displacement_jump[_qp](i)
								+ _zeta * displacement_jump_vel(i)
								- _k1 * -_displacement_jump[_qp](i);
		}
		else
		{
			traction_local(i) = std::pow(_theta ,2) * (1.0 - _damage_old[_qp])
					* _k0 * _displacement_jump[_qp](i)
					+ _zeta * displacement_jump_vel(i);
		}

	}

	Real T_eff;
	T_eff = std::pow(traction_local(0),2);
	T_eff += std::pow(traction_local(1),2) / std::pow(_theta,2);
	T_eff += std::pow(traction_local(2),2) / std::pow(_theta,2);
	_T_eff[_qp] = std::sqrt(T_eff);


	return traction_local;
}


RankTwoTensor
PandolfiOrtizReversible::computeTractionDerivatives()
{
	RankTwoTensor traction_jump_derivatives_local;
	traction_jump_derivatives_local.zero();

	return traction_jump_derivatives_local;
}

