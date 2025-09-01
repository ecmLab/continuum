//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "tecm_testTestApp.h"
#include "tecm_testApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
tecm_testTestApp::validParams()
{
  InputParameters params = tecm_testApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

tecm_testTestApp::tecm_testTestApp(InputParameters parameters) : MooseApp(parameters)
{
  tecm_testTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

tecm_testTestApp::~tecm_testTestApp() {}

void
tecm_testTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  tecm_testApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"tecm_testTestApp"});
    Registry::registerActionsTo(af, {"tecm_testTestApp"});
  }
}

void
tecm_testTestApp::registerApps()
{
  registerApp(tecm_testApp);
  registerApp(tecm_testTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
tecm_testTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  tecm_testTestApp::registerAll(f, af, s);
}
extern "C" void
tecm_testTestApp__registerApps()
{
  tecm_testTestApp::registerApps();
}
