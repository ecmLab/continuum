/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.cc to edit this template
 */

/* 
 * File:   SurfaceConcentration2ParameterAppox.C
 * Author: srinathcs
 * 
 * Created on March 23, 2022, 9:15 AM
 */

#include "SurfaceConcentration2ParameterApprox.h"

registerMooseObject("electro_chemo_mechApp", SurfaceConcentration2ParameterAppox);

InputParameters
SurfaceConcentration2ParameterAppox::validParams()
{
    InputParameters params = ADKernelValue::validParams();
    params.addClassDescription("Residual term (u - v - prop*scale) to set "
           "the variable u equal to the value of couple variable "
            "v - scale * material property");
    params.addParam<std::string>("base_name", "Base name for material class");
    params.addRequiredParam<MaterialPropertyName>(
      "flux_property", "Name of material property to be used in the kernel");
    params.addRequiredParam<MaterialPropertyName>("solid_diffusivity", 
            "Material Property giving solid phase diffusivity");
    params.addRequiredCoupledVar("average_concentration_var", "Variable for averaged concentration");
    params.addRequiredParam<MaterialPropertyName>("particle_size", "Particle size material property");
    params.addParam<Real>("scale",1.0,"SCaling factor");
    return params;
}

SurfaceConcentration2ParameterAppox::SurfaceConcentration2ParameterAppox(const InputParameters & parameters)
        : ADKernelValue(parameters),
        _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
        _coupled_var(adCoupledValue("average_concentration_var")),
        _flux(getADMaterialPropertyByName<Real>(_base_name + getParam<MaterialPropertyName>("flux_property"))),
        _diff(getADMaterialPropertyByName<Real>(getParam<MaterialPropertyName>("solid_diffusivity"))), 
        _particle_size(getADMaterialPropertyByName<Real>(_base_name + getParam<MaterialPropertyName>("particle_size"))), 
        _scale(getParam<Real>("scale"))
        
{  
    
}        

ADReal
SurfaceConcentration2ParameterAppox::precomputeQpResidual()
{
    return -_flux[_qp]/5.0/_diff[_qp] * _particle_size[_qp] * _scale
            + _coupled_var[_qp] - _u[_qp];
}