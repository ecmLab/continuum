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
 * File:   ADChemoMechanoDiffusion.C
 * Author: srinath
 * 
 * Created on October 22, 2020, 11:18 AM
 */

#include "ADChemoMechanoDiffusion.h"


registerMooseObject("electro_chemo_mechApp", ADChemoMechanoDiffusion);

InputParameters
ADChemoMechanoDiffusion::validParams()
{
  auto params = ADChemoMechanoDiffusionTempl<Real>::validParams();
  params.addClassDescription(
      "Diffusion equation kernel that takes an isotropic diffusivity from a material property");

  return params;
}

ADChemoMechanoDiffusion::ADChemoMechanoDiffusion(const InputParameters & parameters)
  : ADChemoMechanoDiffusionTempl<Real>(parameters)
{  
}
