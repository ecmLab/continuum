/*
 * Copyright (C) 2020 srinath
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * File:   ADButlerVolmerBC.h
 * Author: srinath
 * 
 * Created on October 20, 2020, 8:25 AM
 */

#pragma once 

#include "ADIntegratedBC.h"

class ADButlerVolmerBC : public ADIntegratedBC
{
public:
    static InputParameters validParams();
    
    ADButlerVolmerBC(const InputParameters & parameters);
protected:
    ADReal computeQpResidual() override;
    ADReal _i0;
    ADReal _faraday; 
    ADReal _temp; 
    ADReal _gas_constant;
    ADReal _current; 
    const ADMaterialProperty<Real> * _equilibrium_potential; 
    const Function * _func;

};