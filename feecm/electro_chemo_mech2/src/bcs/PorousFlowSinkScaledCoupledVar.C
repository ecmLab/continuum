/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   PorousFlowSinkScaledCoupledVar.C
 * Author: srinath
 * 
 * Created on October 19, 2021, 6:36 PM
 */

#include "PorousFlowSinkScaledCoupledVar.h"

registerMooseObject("electro_chemo_mechApp", PorousFlowSinkScaledCoupledVar);

InputParameters
PorousFlowSinkScaledCoupledVar::validParams()
{
    InputParameters params = PorousFlowSink::validParams();
    params.addRequiredCoupledVar("v", "Coupled variable setting the gradient on the boundary.");
    params.addParam<Real>("scale", "Scaling factor for variable to impose BC");
    params.addClassDescription("Imposes the integrated boundary condition "
                             "$\\frac{\\partial u}{\\partial n}=v$, "
                             "where $v$ is a variable.");
    return params;
}

PorousFlowSinkScaledCoupledVar::PorousFlowSinkScaledCoupledVar(const InputParameters& parameters) 
    : PorousFlowSink(parameters),
    _scale(isParamValid("scale") ? getParam<Real>("scale") : 1.0),
    _coupled_var(coupledValue("v"))
     
{     
}

Real
PorousFlowSinkScaledCoupledVar::computeQpResidual()
{
    return PorousFlowSink::computeQpResidual() * _coupled_var[_qp] * _scale;
}
