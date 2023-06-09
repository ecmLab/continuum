//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "BoundaryBase.h"

/**
 * Base class for heat structure boundary components
 */
class HSBoundary : public BoundaryBase
{
public:
  HSBoundary(const InputParameters & params);

  virtual void check() const override;

protected:
  /// Boundary names for which the boundary component applies
  const std::vector<BoundaryName> & _boundary;

  /// Heat structure name
  const std::string & _hs_name;

public:
  static InputParameters validParams();
};
