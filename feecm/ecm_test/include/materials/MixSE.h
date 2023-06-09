
#pragma once

#include "ADMaterial.h"
class MixSE;

/**
 * Material-derived objects override the computeQpProperties()
 * function.  They must declare and compute material properties for
 * use by other objects in the calculation such as Kernels and
 * BoundaryConditions.
 */
class MixSE : public ADMaterial
{
public:
  static InputParameters validParams();

  MixSE(const InputParameters & parameters);

protected:
  /// Necessary override. This is where the values of the properties are computed.
  virtual void computeQpProperties() override;

  /// Values from the input file
  const Real & _inIonicConductivity;
  const Real & _inMetalConductivity;
  const Real & _inElectronicConductivity;
  const Real & _inInletCurrent;
  const Real & _inExchangeCurrent;
  const Real & _inReactionRate;
  const Real & _inElectronConcentration;
  const Real & _inLiPotAnode;
  const Real & _inLiPotCathode;

  /// The ionic conductivity of Li+ in SE
  ADMaterialProperty<Real> & _ionic_conductivity;
  
  /// The electric conductivity of electron in Li metal
  ADMaterialProperty<Real> & _metal_conductivity;
  
  /// The electronic conductivity of SE
  ADMaterialProperty<Real> & _electronic_conductivity;

  /// The inlet current density from cathode boundary
  ADMaterialProperty<Real> & _inlet_current;

  /// The exchange current density at Li metal/SE interface
  ADMaterialProperty<Real> & _exchange_current;
  
  /// The reaction rate of the Li/Li+ reaction
  ADMaterialProperty<Real> & _reaction_rate;

  /// The relative electron concentration in the SE
  ADMaterialProperty<Real> & _electron_concentration;

  /// The Li chemical potential in Anode, unit V.
  ADMaterialProperty<Real> & _LiPotAnode;

  /// The Li chemical potential in Cathode, unit V.
  ADMaterialProperty<Real> & _LiPotCathode;

//  usingMaterialMembers;
//  using ADMaterial<compute_stage>::_communicator;
};
