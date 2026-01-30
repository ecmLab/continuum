//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "mece789TestApp.h"
#include "mece789App.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
mece789TestApp::validParams()
{
  InputParameters params = mece789App::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

mece789TestApp::mece789TestApp(InputParameters parameters) : MooseApp(parameters)
{
  mece789TestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

mece789TestApp::~mece789TestApp() {}

void
mece789TestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  mece789App::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"mece789TestApp"});
    Registry::registerActionsTo(af, {"mece789TestApp"});
  }
}

void
mece789TestApp::registerApps()
{
  registerApp(mece789App);
  registerApp(mece789TestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
mece789TestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  mece789TestApp::registerAll(f, af, s);
}
extern "C" void
mece789TestApp__registerApps()
{
  mece789TestApp::registerApps();
}
