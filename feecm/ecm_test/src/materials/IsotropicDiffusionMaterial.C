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
 * File:   IsotropicDiffusionMaterial.C
 * Author: srinath
 * 
 * Created on October 22, 2020, 9:53 AM
 */

#include "IsotropicDiffusionMaterial.h"

registerADMooseObject("electro_chemo_mechApp", ADIsotropicDiffusionMaterial);

InputParameters ADIsotropicDiffusionMaterial::validParams()
{
    InputParameters params = Material::validParams();
    params.addClassDescription("Provides a diffusivity tensor based on 3 principal directions");
    params.addParam<std::string>("base_name", "Base name for material class");
    params.addRequiredParam<Real>("diffusion_coef", "Diffusivity");
    params.addParam<MaterialPropertyName>("diffusivity_name", "diffusivity", 
            "Name of the diffusivity");
    return params;
}

ADIsotropicDiffusionMaterial::ADIsotropicDiffusionMaterial
    (const InputParameters & parameters) : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _diffusion_coef(getParam<Real>("diffusion_coef")),
    _diffusivity(declareADProperty<Real>(_base_name + getParam<MaterialPropertyName>("diffusivity_name")))
{
    
}

void ADIsotropicDiffusionMaterial::computeQpProperties()
{
    
    _diffusivity[_qp] = computeQpCorrection() * _diffusion_coef;
}
