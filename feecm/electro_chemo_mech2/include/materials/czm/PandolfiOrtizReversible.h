/*
 * PandolfiOrtizReversible.h
 *
 *  Created on: Jun 17, 2020
 *      Author: srinath
 */

#pragma once

#include "CZMMaterialBase.h"

/**
 * Implementation of Pandolfi and Ortiz non stateful traction-separation law
 * Ortiz and Pandolfi, 1999
 * The traction separation law is linear in loading and unloading
 */

class PandolfiOrtizReversible : public CZMMaterialBase {
public:
	static InputParameters validParams();
	PandolfiOrtizReversible(const InputParameters & parameters);

protected:
	virtual void initQpStatefulProperties() override;

	virtual RealVectorValue computeTraction() override;

    virtual RankTwoTensor computeTractionDerivatives() override;



  /// interface displacement jump when traction reaches max stress
  const Real _d1;

  /// interface displacement jump when traction reaches zero
  const Real _d2;

  /// Compressive stiffness in normal opening mode
  const Real _k1;

  /// Tensile stiffness of undamaged interface
  const Real _k0;

  /// Relative stiffness of interface between shear and tension
  const Real _theta;

  /// Artficial viscosity for numerical stabilization
  const Real _zeta;

  /// Max. initial interface strength
  const Real _sigma_max;

  /// Effective separation and Traction
  MaterialProperty<Real> & _delta_eff;
  const MaterialProperty<Real> & _delta_eff_old;
  MaterialProperty<Real> & _T_eff;
  const MaterialProperty<Real> & _T_eff_old;

  /// Damage
  MaterialProperty<Real> & _damage;
  const MaterialProperty<Real> & _damage_old;

  /// Damage modified interface strength
  MaterialProperty<Real> & _damage_strength;
  const MaterialProperty<Real> & _damage_strength_old;



};
