/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.cc to edit this template
 */

/* 
 * File:   ADMultipleButlerVolmerForce.C
 * Author: srinath
 * 
 * Created on March 24, 2022, 4:48 PM
 */

#include <vector>

#include "ADMultipleButlerVolmerForce.h"

registerADMooseObject("electro_chemo_mechApp", ADMultipleButlerVolmerForce);

InputParameters
ADMultipleButlerVolmerForce::validParams()
{
    InputParameters params = ADKernel::validParams();
    params.addClassDescription("Implements a body force that is dependent on a material property");
    params.addParam<std::string>("base_name", "Base Name");
      params.addRequiredParam<std::vector<MaterialName>>(
      "materials",
      "The material objects to use to calculate kinetic source term");
    params.addParam<Real>("scale", "Scaling coefficient");
    params.addParam<MaterialPropertyName>("mat_prop_name", "Name of material property");
    return params;
}

ADMultipleButlerVolmerForce::ADMultipleButlerVolmerForce(const InputParameters & parameters)
        : ADKernel(parameters),
        _scale(getParam<Real>("scale")),
        _num_materials(getParam<std::vector<MaterialName>>("materials").size())
{
    std::vector<MaterialName> models = getParam<std::vector<MaterialName>>("materials");
    for (unsigned int i = 0; i < _num_materials; ++i) {
        ADButlerVolmerMaterial * xxx = dynamic_cast<ADButlerVolmerMaterial *>(&this->getMaterialByName(models[i], true));
        if (xxx)
        {
//            _materials.push_back(xxx);
            auto base_name = (xxx->isParamValid("base_name") ? xxx->getParam<std::string>("base_name") + "_" : "");
            auto abc = base_name + getParam<MaterialPropertyName>("mat_prop_name");
            auto xyz = xxx->hasADMaterialPropertyByName<Real>(abc);
            if (!xyz)
                mooseError("Material " + models[i] + " requires ", 
                        getParam<MaterialPropertyName>("mat_prop_name"));
            else{
                _mat_properties.push_back(&(xxx->getADMaterialProperty<Real>(abc)));
            }
        } else
        {
            mooseError("Material " + models[i] + 
                    " is not compatible with ADMultiCoupledForce");
        }
    }
    
}

ADReal
ADMultipleButlerVolmerForce::computeQpResidual()
{
    ADReal res = 0.0;
    for (unsigned int i = 0; i < _num_materials; ++i) {
//        auto m = _mat_properties[i];
        res += _scale * (*_mat_properties[i])[_qp];
    }
    return  res * _test[_i][_qp];
}
