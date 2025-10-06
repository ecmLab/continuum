//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "edl_solidTestApp.h"
#include "edl_solidApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
edl_solidTestApp::validParams()
{
  InputParameters params = edl_solidApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

edl_solidTestApp::edl_solidTestApp(InputParameters parameters) : MooseApp(parameters)
{
  edl_solidTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

edl_solidTestApp::~edl_solidTestApp() {}

void
edl_solidTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  edl_solidApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"edl_solidTestApp"});
    Registry::registerActionsTo(af, {"edl_solidTestApp"});
  }
}

void
edl_solidTestApp::registerApps()
{
  registerApp(edl_solidApp);
  registerApp(edl_solidTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
edl_solidTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  edl_solidTestApp::registerAll(f, af, s);
}
extern "C" void
edl_solidTestApp__registerApps()
{
  edl_solidTestApp::registerApps();
}
