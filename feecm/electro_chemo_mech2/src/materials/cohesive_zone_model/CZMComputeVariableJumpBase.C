/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CZMComputeVariableJumpBase.C
 * Author: srinath
 * 
 * Created on January 25, 2022, 1:06 PM
 */

#include "CZMComputeVariableJumpBase.h"
#include "CohesiveZoneModelTools.h"

InputParameters
CZMComputeVariableJumpBase::validParams()
{
  InputParameters params = InterfaceMaterial::validParams();
  params.addClassDescription("Base class used to compute the variable jump across a czm "
                             "interface in local coordinates");
  params.addRequiredCoupledVar("variable",
                               "The string of variable suitable for the problem statement");
  params.addCoupledVar("neighbor_variable", "Optinal name of neighbor variable");
  params.suppressParameter<bool>("use_displaced_mesh");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addParam<bool>("include_gap", true, "Gap inclusion");
 
  return params;
}

CZMComputeVariableJumpBase::CZMComputeVariableJumpBase(const InputParameters & parameters)
  : InterfaceMaterial(parameters),
    _base_name(isParamValid("base_name") && !getParam<std::string>("base_name").empty()
                   ? getParam<std::string>("base_name") + "_"
                   : ""),
    _include_gap(getParam<bool>("include_gap")),
     _czm_total_rotation(_include_gap ? 
         &getMaterialPropertyByName<RankTwoTensor>(_base_name + "czm_total_rotation") : nullptr),
     _total_rotation(
        declarePropertyByName<RankTwoTensor>(_base_name + "total_rotation")),
     _variable_jump_global(
        declarePropertyByName<Real>(_base_name + "variable_jump_global")),
    _interface_variable_jump(
        declarePropertyByName<Real>(_base_name + "interface_variable_jump"))

{
    _variable = &coupledValue("variable",0);
    if (isParamValid("neighbor_variable"))
        _variable_neighbor = &coupledNeighborValue("neighbor_variable",0);
    else 
        _variable_neighbor = &coupledNeighborValue("variable",0);
}

void
CZMComputeVariableJumpBase::computeQpProperties()
{
    if (_include_gap && _czm_total_rotation)
        _total_rotation[_qp] = (*_czm_total_rotation)[_qp];
    else{
        _total_rotation[_qp] = CohesiveZoneModelTools::computeReferenceRotation(_normals[_qp], _mesh.dimension());
    }
  // computing the displacement jump
    _variable_jump_global[_qp] = (*_variable_neighbor)[_qp] - (*_variable)[_qp];

  computeLocalVariableJump();
}