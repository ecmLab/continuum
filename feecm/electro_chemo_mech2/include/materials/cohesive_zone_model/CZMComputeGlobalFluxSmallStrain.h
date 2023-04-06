/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CZMComputeGlobalFlux.h
 * Author: srinath
 *
 * Created on January 25, 2022, 2:27 PM
 */



#pragma once

#include "CZMcomputeGlobalFluxBase.h"
/**
 * This class uses the interface traction and its derivatives w.r.t. the interface displacment jump
 * to compute their respective values in global coordinates under the small strain assumption. The
 * values computed by this object are used by the CZMInterfaceKernelSmallStrain to add the proper
 * residual to the system and to compute the analytic Jacobian.
 */

class CZMComputeGlobalFluxSmallStrain : public CZMcomputeGlobalFluxBase
{
public:
  static InputParameters validParams();
  CZMComputeGlobalFluxSmallStrain(const InputParameters & parameters);

protected:
  /// method used to compute the traction and it's derivatives in global coordinates.
  void computeEquilibriumFluxAndDerivatives() override;
};
