//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "electro_chemo_mech2TestApp.h"
#include "electro_chemo_mech2App.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
electro_chemo_mech2TestApp::validParams()
{
  InputParameters params = electro_chemo_mechApp::validParams();
  return params;
}

electro_chemo_mech2TestApp::electro_chemo_mech2TestApp(InputParameters parameters) : MooseApp(parameters)
{
  electro_chemo_mech2TestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

electro_chemo_mech2TestApp::~electro_chemo_mech2TestApp() {}

void
electro_chemo_mech2TestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  electro_chemo_mechApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"electro_chemo_mech2TestApp"});
    Registry::registerActionsTo(af, {"electro_chemo_mech2TestApp"});
  }
}

void
electro_chemo_mech2TestApp::registerApps()
{
  registerApp(electro_chemo_mechApp);
  registerApp(electro_chemo_mech2TestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
electro_chemo_mech2TestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  electro_chemo_mech2TestApp::registerAll(f, af, s);
}
extern "C" void
electro_chemo_mech2TestApp__registerApps()
{
  electro_chemo_mech2TestApp::registerApps();
}
