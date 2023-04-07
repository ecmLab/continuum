//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "ecmMechanicsTestApp.h"
#include "ecmMechanicsApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
ecmMechanicsTestApp::validParams()
{
  InputParameters params = ecmMechanicsApp::validParams();
  return params;
}

ecmMechanicsTestApp::ecmMechanicsTestApp(InputParameters parameters) : MooseApp(parameters)
{
  ecmMechanicsTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

ecmMechanicsTestApp::~ecmMechanicsTestApp() {}

void
ecmMechanicsTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  ecmMechanicsApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"ecmMechanicsTestApp"});
    Registry::registerActionsTo(af, {"ecmMechanicsTestApp"});
  }
}

void
ecmMechanicsTestApp::registerApps()
{
  registerApp(ecmMechanicsApp);
  registerApp(ecmMechanicsTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
ecmMechanicsTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ecmMechanicsTestApp::registerAll(f, af, s);
}
extern "C" void
ecmMechanicsTestApp__registerApps()
{
  ecmMechanicsTestApp::registerApps();
}
