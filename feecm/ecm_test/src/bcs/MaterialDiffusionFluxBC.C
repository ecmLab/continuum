/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MaterialDiffusionFluxBC.C
 * Author: srinath
 * 
 * Created on February 12, 2021, 5:22 PM
 */

#include "MaterialDiffusionFluxBC.h"


registerMooseObject("electro_chemo_mechApp", MaterialDiffusionFluxBC);

InputParameters
MaterialDiffusionFluxBC::validParams()
{
    InputParameters params = ADIntegratedBC::validParams();
    params.addParam<MaterialPropertyName>("diffusivity", "diffusivity", "Name of diffusivity");
    return params;
}


MaterialDiffusionFluxBC::MaterialDiffusionFluxBC(const InputParameters& parameters)
        : ADIntegratedBC(parameters),
        _diffusivity(&getADMaterialProperty<Real>("diffusivity"))
{
    
}


ADReal
MaterialDiffusionFluxBC::computeQpResidual()
{
    return -(*_diffusivity)[_qp] * _grad_u[_qp] * _normals[_qp] * _test[_i][_qp];
}