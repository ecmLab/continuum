/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CZMComputeGlobalFlux.C
 * Author: srinath
 * 
 * Created on January 25, 2022, 2:27 PM
 */

#include "CZMComputeGlobalFluxSmallStrain.h"


registerMooseObject("electro_chemo_mechApp", CZMComputeGlobalFluxSmallStrain);

InputParameters
CZMComputeGlobalFluxSmallStrain::validParams()
{
  InputParameters params = CZMcomputeGlobalFluxBase::validParams();

  params.addClassDescription(
      "Computes the czm flux in global coordinates for a small strain kinematic formulation");
  return params;
}

CZMComputeGlobalFluxSmallStrain::CZMComputeGlobalFluxSmallStrain(
    const InputParameters & parameters)
  : CZMcomputeGlobalFluxBase(parameters)
{
}

void
CZMComputeGlobalFluxSmallStrain::computeEquilibriumFluxAndDerivatives()
{
  _flux_global[_qp] = _interface_flux[_qp];
  _dflux_dvariablejump_global[_qp] = _dinterface_flux_dvariablejump[_qp];
  _dflux_djump_global[_qp] = _total_rotation[_qp] * _dinterface_flux_djump[_qp];
                              
}
