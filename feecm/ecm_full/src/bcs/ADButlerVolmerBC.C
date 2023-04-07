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
 * File:   ADButlerVolmerBC.C
 * Author: srinath
 * 
 * Created on October 20, 2020, 8:25 AM
 */

#include "ADButlerVolmerBC.h"

#include "Function.h"


registerMooseObject("electro_chemo_mechApp", ADButlerVolmerBC);

InputParameters
ADButlerVolmerBC::validParams()
{
    InputParameters params = ADIntegratedBC::validParams();
    params.addRequiredParam<Real>("exchange_current_density", "Value of exchange"
            "current density");
    params.addParam<Real>("faraday", 96485.3329, "Faraday's Constant");
    params.addParam<Real>("R", 8.3145, "Universal Gas Constant");
    params.addParam<Real>("Temperature", 298, "Value of temperature to use");
    params.addParam<Real>("current_density", "Applied Current Density");
    params.addParam<FunctionName>("current_density_function", "", "Function for current density");
    return params;
}

ADButlerVolmerBC::ADButlerVolmerBC(const InputParameters & parameters)
        : ADIntegratedBC(parameters),
        _i0(getParam<Real>("exchange_current_density")),
        _faraday(getParam<Real>("faraday")),
        _gas_constant(getParam<Real>("R")),
        _temp(getParam<Real>("Temperature")),
        _current(isParamValid("current_density") ? getParam<Real>("current_density") : 0), 
        _equilibrium_potential(&getADMaterialProperty<Real>("equilibrium_potential") ? &getADMaterialProperty<Real>("equilibrium_potential") : NULL), 
        _func(getParam<FunctionName>("current_density_function") != "" 
                            ? &getFunction("current_density_function") : NULL)
{
    if (isParamValid("current_density") && _func)
    {
        mooseError("Cannot define both current density and current density function");
    }
}

ADReal
ADButlerVolmerBC::computeQpResidual()
{
    auto prefac = 0.5*_faraday/(_gas_constant * _temp);
    auto deltaPhi = prefac * _u[_qp];
    if (_equilibrium_potential)
        deltaPhi -= prefac * (*_equilibrium_potential)[_qp];
    auto fac = std::exp(deltaPhi) - std::exp(-deltaPhi);
    ADReal i = 0;
    if (isParamValid("current_density"))
        i = _current;
    if (_func)
        i = _func->value(_t, _q_point[_qp]);
    return (i - _i0 * fac) * _test[_i][_qp];
}