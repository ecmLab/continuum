//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "ecmElectrochemTestApp.h"
#include "ecmElectrochemApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
ecmElectrochemTestApp::validParams()
{
  InputParameters params = ecmElectrochemApp::validParams();
  return params;
}

ecmElectrochemTestApp::ecmElectrochemTestApp(InputParameters parameters) : MooseApp(parameters)
{
  ecmElectrochemTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

ecmElectrochemTestApp::~ecmElectrochemTestApp() {}

void
ecmElectrochemTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  ecmElectrochemApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"ecmElectrochemTestApp"});
    Registry::registerActionsTo(af, {"ecmElectrochemTestApp"});
  }
}

void
ecmElectrochemTestApp::registerApps()
{
  registerApp(ecmElectrochemApp);
  registerApp(ecmElectrochemTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
ecmElectrochemTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ecmElectrochemTestApp::registerAll(f, af, s);
}
extern "C" void
ecmElectrochemTestApp__registerApps()
{
  ecmElectrochemTestApp::registerApps();
}
