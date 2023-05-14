//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "ecBetaTestApp.h"
#include "ecBetaApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
ecBetaTestApp::validParams()
{
  InputParameters params = ecBetaApp::validParams();
  return params;
}

ecBetaTestApp::ecBetaTestApp(InputParameters parameters) : MooseApp(parameters)
{
  ecBetaTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

ecBetaTestApp::~ecBetaTestApp() {}

void
ecBetaTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  ecBetaApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"ecBetaTestApp"});
    Registry::registerActionsTo(af, {"ecBetaTestApp"});
  }
}

void
ecBetaTestApp::registerApps()
{
  registerApp(ecBetaApp);
  registerApp(ecBetaTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
ecBetaTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ecBetaTestApp::registerAll(f, af, s);
}
extern "C" void
ecBetaTestApp__registerApps()
{
  ecBetaTestApp::registerApps();
}
