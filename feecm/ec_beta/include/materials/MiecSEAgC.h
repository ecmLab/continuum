
#pragma once

#include "ADMaterial.h"
/**
 * Material-derived objects override the computeQpProperties()
 * function.  They must declare and compute material properties for
 * use by other objects in the calculation such as Kernels and
 * BoundaryConditions.
 */
class MiecSEAgC : public ADMaterial
{
public:
  static InputParameters validParams();

  MiecSEAgC(const InputParameters & parameters);

protected:
  /// Necessary override. This is where the values of the properties are computed.
  virtual void computeQpProperties() override;

  /// Values from the input file
  const Real & _inIonicConductivitySE;
  const Real & _inIonicConductivityAgC;
  const Real & _inElectronicConductivitySE;
  const Real & _inElectronicConductivityAgC;
  const Real & _inInletCurrent;
  const Real & _inExchangeCurrent;
  const Real & _inReactionRate;
  const Real & _inElectronConcentration;
  const Real & _inLiPotAnode;
  const Real & _inLiPotCathode;

  /// The ionic conductivity of Li+ in SE
  ADMaterialProperty<Real> & _ionic_conductivity_SE;
  ADMaterialProperty<Real> & _ionic_conductivity_AgC;
  
  /// The electronic conductivity of SE
  ADMaterialProperty<Real> & _electronic_conductivity_SE;
  ADMaterialProperty<Real> & _electronic_conductivity_AgC;

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
