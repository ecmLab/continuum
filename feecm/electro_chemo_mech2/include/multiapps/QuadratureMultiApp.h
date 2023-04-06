/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   QuadratureMultiApp.h
 * Author: srinathcs
 *
 * Created on March 22, 2022, 2:34 PM
 */

#pragma once

#include "TransientMultiApp.h"
#include "BlockRestrictable.h"

/**
 * Automatically generates Sub-App positions from centroids of elements in the master mesh.
 */
class QuadratureMultiApp : public TransientMultiApp, public BlockRestrictable
{
public:
  static InputParameters validParams();

  QuadratureMultiApp(const InputParameters & parameters);

protected:
  /**
   * fill in _positions with the positions of the sub-aps
   */
  virtual void fillPositions() override;
};
