/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ADRankTwoTensorNormalToBoundary.C
 * Author: srinath
 * 
 * Created on January 13, 2021, 12:46 PM
 */

#include "ADRankTwoTensorNormalToBoundary.h"
#include "Assembly.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("electro_chemo_mechApp", ADRankTwoTensorNormalToBoundary);

// defineLegacyParams(ADRankTwoTensorNormalToBoundary);

InputParameters
ADRankTwoTensorNormalToBoundary :: validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Compute scalar value of rank two tensor normal to a boundary");
  params.addRequiredParam<MaterialPropertyName>(
      "tensor",
      "The name of the RankTwoTensor Material Property");
  return params;
    
}        

ADRankTwoTensorNormalToBoundary::ADRankTwoTensorNormalToBoundary(const InputParameters & parameters)
        : AuxKernel(parameters),
        _normals(_assembly.normals()),
        _tensor(getADMaterialProperty<RankTwoTensor>("tensor"))
{    
}
 
Real 
ADRankTwoTensorNormalToBoundary::computeValue()
{
    auto temp1 = _normals[_qp] * _tensor[_qp].column(0);
    auto temp2 = _normals[_qp] * _tensor[_qp].column(1);
    auto temp3 = _normals[_qp] * _tensor[_qp].column(2);
    Point p1 = (MetaPhysicL::raw_value(temp1), MetaPhysicL::raw_value(temp2), MetaPhysicL::raw_value(temp3));
    Point p2 = (MetaPhysicL::raw_value(_normals[_qp](0)), MetaPhysicL::raw_value(_normals[_qp](1)), MetaPhysicL::raw_value(_normals[_qp](2)));
    return p1 * p2;
//    return MetaPhysicL::raw_value(_normals[_qp].transpose() * _tensor[_qp] * _normals[_qp]);

}

