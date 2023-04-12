//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "ecm_fullTestApp.h"
#include "ecm_fullApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
ecm_fullTestApp::validParams()
{
  InputParameters params = electro_chemo_mechApp::validParams();
  return params;
}

ecm_fullTestApp::ecm_fullTestApp(InputParameters parameters) : MooseApp(parameters)
{
  ecm_fullTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

ecm_fullTestApp::~ecm_fullTestApp() {}

void
ecm_fullTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  electro_chemo_mechApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"ecm_fullTestApp"});
    Registry::registerActionsTo(af, {"ecm_fullTestApp"});
  }
}

void
ecm_fullTestApp::registerApps()
{
  registerApp(ecm_fullApp);
  registerApp(ecm_fullTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
ecm_fullTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ecm_fullTestApp::registerAll(f, af, s);
}
extern "C" void
ecm_fullTestApp__registerApps()
{
  cm_fullTestApp::registerApps();
}
