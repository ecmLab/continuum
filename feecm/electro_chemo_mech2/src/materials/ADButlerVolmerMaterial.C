/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.cc to edit this template
 */

/* 
 * File:   ADButlerVolmerMaterial.C
 * Author: srinath
 * 
 * Created on March 18, 2022, 6:41 AM
 */

#include "ADButlerVolmerMaterial.h"

registerADMooseObject("electro_chemo_mechApp", ADButlerVolmerMaterial);

InputParameters
ADButlerVolmerMaterial::validParams()
{
    InputParameters params = ADMaterial::validParams();
    params.addClassDescription("Base class for all Butler-Volmer Materials");
    params.addParam<std::string>("base_name", "Base name for material class");
    params.addParam<Real>("faraday", 96.4853329, "Faraday's Constant");
    params.addParam<Real>("gas_constant", 8.3145, "Universal Gas Constant");
    params.addParam<Real>("temperature", 298, "Value of temperature to use");
    params.addParam<bool>("include_equilibrium", false, 
            "Choose whether to include equilibrium_potential");  
    params.addRequiredCoupledVar("electrolyte_potential_var", 
            "Name of Variable for the electronic potential");
    params.addRequiredCoupledVar("electrode_potential_var", 
            "Name of Variable for the ionic potential");
    params.addRequiredParam<MaterialPropertyName>("exchange_current_density", "Exchange Current Density Material Name");
    params.addRequiredParam<MaterialPropertyName>("surface_to_volume", "Surface to volume ratio material property");
    params.addParam<MaterialPropertyName>("volume_fraction", "volume_fraction", "Volume fraction of particle");
    params.addParam<Real>("specific_area", 1.0, "Surface to volume ratio for the Electrode");
    return params;
}


ADButlerVolmerMaterial::ADButlerVolmerMaterial(const InputParameters & parameters)
        : ADMaterial(parameters),
        _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
        _gas_constant(getParam<Real>("gas_constant")),
        _temperature(getParam<Real>("temperature")),
        _faraday(getParam<Real>("faraday")),
        _electrolyte_potential(isParamValid("electrolyte_potential_var") ?  
                               &adCoupledValue("electrolyte_potential_var") : 
                                nullptr), 
        _electrode_potential(isParamValid("electrode_potential_var") ?  
                               &adCoupledValue("electrode_potential_var") : 
                                nullptr),
        _include_equil(getParam<bool>("include_equilibrium")), 
        _equilibirium_potential(_include_equil 
                            ? &getADMaterialPropertyByName<Real>(_base_name + "equilibrium_potential")
                            : nullptr), 
        _bvcurrent(declareADPropertyByName<Real>(_base_name + "butler_volmer_current")),
        _bvflux(declareADPropertyByName<Real>(_base_name + "butler_volmer_flux")),
        _bvflux_force(declareADPropertyByName<Real>(_base_name + "butler_volmer_flux_force")),
        _bvcurrent_force(declareADProperty<Real>(_base_name + "butler_volmer_current_force")),
        _eta(declareADProperty<Real>(_base_name + "local_overpotential")),
        _surface_to_volume(getADMaterialPropertyByName<Real>(_base_name +  getParam<MaterialPropertyName>("surface_to_volume"))), 
        _volume_fraction(getADMaterialPropertyByName<Real>(_base_name + getParam<MaterialPropertyName>("volume_fraction"))),
        _i0(&getADMaterialPropertyByName<Real>(_base_name +  getParam<MaterialPropertyName>("exchange_current_density"))), 
        _a0(getParam<Real>("specific_area"))
{
    
}

void
ADButlerVolmerMaterial::computeQpProperties()
{
    auto RTF = _gas_constant * _temperature / _faraday;
    auto eta = (*_electrode_potential)[_qp] - (*_electrolyte_potential)[_qp];
    if ( _equilibirium_potential) eta -= (*_equilibirium_potential)[_qp];
    _eta[_qp] = eta;               
    _bvcurrent[_qp] = (*_i0)[_qp] * 2.0 * std::sinh(0.5*eta/RTF);
    _bvflux[_qp] = _bvcurrent[_qp]/_faraday;
    _bvflux_force[_qp] = _bvflux[_qp] * _surface_to_volume[_qp];
    _bvcurrent_force[_qp] = _bvcurrent[_qp] * _surface_to_volume[_qp] * _volume_fraction[_qp];
}
