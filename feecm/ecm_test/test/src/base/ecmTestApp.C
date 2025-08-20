//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "ecmTestApp.h"
#include "ecmApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
ecmTestApp::validParams()
{
  InputParameters params = ecmApp::validParams();
  return params;
}

ecmTestApp::ecmTestApp(InputParameters parameters) : MooseApp(parameters)
{
  ecmTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

ecmTestApp::~ecmTestApp() {}

void
ecmTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  ecmApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"ecmTestApp"});
    Registry::registerActionsTo(af, {"ecmTestApp"});
  }
}

void
ecmTestApp::registerApps()
{
  registerApp(ecmApp);
  registerApp(ecmTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
ecmTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ecmTestApp::registerAll(f, af, s);
}
extern "C" void
ecmTestApp__registerApps()
{
  ecmTestApp::registerApps();
}
