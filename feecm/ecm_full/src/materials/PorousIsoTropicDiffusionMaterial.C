/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.cc to edit this template
 */

/* 
 * File:   PorousIsoTropicDiffusionMaterial.C
 * Author: srinathcs
 * 
 * Created on March 23, 2022, 4:25 PM
 */

#include "PorousIsoTropicDiffusionMaterial.h"

registerMooseObject("electro_chemo_mechApp", PorousIsotropicDiffusionMaterial);

InputParameters
PorousIsotropicDiffusionMaterial::validParams()
{
    InputParameters params = ADIsotropicDiffusionMaterial::validParams();
    params.addClassDescription("Compute the Effective diffusive "
            "or conductive properties of porous medium");
    params.addRequiredParam<Real>("volume_fraction", "Volume fraction of material");
    params.addParam<Real>("Bruggeman_factor", 1.5, "Exponenent of Bruggeman equation");
    MooseEnum diffusionModel("BRUGGEMAN TORTUOSITY USER_DEFINED", "BRUGGEMAN");
    params.addParam<MooseEnum>("diffusion_model", diffusionModel, "Type of correction");
    params.addParam<Real>("tortuosity", 1.0, "Tortuosity of material");
    params.addParam<Real>("correction_factor", 1.0, "Correction factor for user defined material");
    return params;
}

PorousIsotropicDiffusionMaterial::PorousIsotropicDiffusionMaterial(const InputParameters & parameters)
        : ADIsotropicDiffusionMaterial(parameters),
        _volume_fraction(getParam<Real>("volume_fraction")), 
        _bruggeman_factor(getParam<Real>("Bruggeman_factor")),
        _diffusion_model(getParam<MooseEnum>("diffusion_model").getEnum<DiffusionModel>()),
        _tortuosity(getParam<Real>("tortuosity")),
        _correction_factor(getParam<Real>("correction_factor"))
{
    if (_diffusion_model == DiffusionModel::Tortuosity)
    {
        if (!isParamValid("tortuosity"))
            mooseError("Cannot have a tortuosity correction without a tortuosity value");
    }
 
}

Real
PorousIsotropicDiffusionMaterial::computeQpCorrection()
{
    switch(_diffusion_model)
    {
        case DiffusionModel::Bruggeman:
            return std::pow(_volume_fraction, _bruggeman_factor);
        case DiffusionModel::Tortuosity:
            return _volume_fraction/_tortuosity;
        case DiffusionModel::User_Defined:
            return _correction_factor;
    }
}