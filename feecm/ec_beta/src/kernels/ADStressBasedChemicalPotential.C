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
 * File:   ADStressBasedChemicalPotential.C
 * Author: srinath
 * 
 * Created on October 19, 2020, 1:22 PM
 */

#include "ADStressBasedChemicalPotential.h"

registerMooseObject("ecBetaApp", ADStressBasedChemicalPotential);

InputParameters
ADStressBasedChemicalPotential::validParams()
{
    InputParameters params = ADKernel::validParams();
    params.addClassDescription("Stress based chemical potential, calculates the "
            "the contribution of elastic energy and swelling based energy to the "
            "chemical potential. These quantities are all calculated in the reference"
            "space");
    params.addParam<std::string>("base_name", "Material property base name");
    params.set<bool>("use_displaced_mesh", false);
    return params;
}

ADStressBasedChemicalPotential::ADStressBasedChemicalPotential(const InputParameters & parameters)
    : ADKernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _swelling_chemical_potential(getADMaterialProperty<Real>(_base_name + "swelling_chemical_potential")),
    _elastic_chemical_potential(getADMaterialProperty<Real>(_base_name + "elastic_chemical_potential"))
{
    
}

ADReal
ADStressBasedChemicalPotential::computeQpResidual()
{
    return (_u[_qp] - _swelling_chemical_potential[_qp] - 
            _elastic_chemical_potential[_qp]) * _test[_i][_qp];
}