/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   DiffusionAlongPrincipalDirectionsMaterial.C
 * Author: srinath
 * 
 * Created on August 14, 2020, 12:43 PM
 */

#include "DiffusionAlongPrincipalDirectionsMaterial.h"

registerADMooseObject("electro_chemo_mechApp", ADDiffusionAlongPrincipalDirectionsMaterial);

InputParameters ADDiffusionAlongPrincipalDirectionsMaterial::validParams()
{
    InputParameters params = Material::validParams();
    params.addParam<std::string>("base_name", "Base name for material class");
    params.addClassDescription("Provides a diffusivity tensor based on 3 principal directions");
    params.addRequiredParam<RealVectorValue>("diffusivity_vector", "Diffusivity along 3 principal directions ");
    params.addParam<MaterialPropertyName>("diffusivity", "diffusivity", 
            "Name of the diffusivity");
    return params;
}

ADDiffusionAlongPrincipalDirectionsMaterial::ADDiffusionAlongPrincipalDirectionsMaterial
    (const InputParameters & parameters) : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _diffusivity_values(getParam<RealVectorValue>("diffusivity_vector")),
    _swell_normal(getADMaterialProperty<RealVectorValue>("swell_normal")),
    _diffusivity(declareADPropertyByName<RealTensorValue>(_base_name + getParam<MaterialPropertyName>("diffusivity")))
{
    
}

void ADDiffusionAlongPrincipalDirectionsMaterial::computeQpProperties()
{
    ADRealVectorValue m_v1, m_v2, m_v3;
    m_v1 = _swell_normal[_qp];
    m_v3 = {0, 0, 1};
    m_v2 = m_v1.cross(m_v3);
    ADRankTwoTensor D1, D2, D3;

    D1.vectorOuterProduct(m_v1, m_v1);
    D2.vectorOuterProduct(m_v2, m_v2);
    D3.vectorOuterProduct(m_v3, m_v3);  
    
    _diffusivity[_qp] = _diffusivity_values(0)* D1 + 
                        _diffusivity_values(1)* D2 +
                        _diffusivity_values(2)* D3;
}
