#include "mece789App.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
mece789App::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

mece789App::mece789App(InputParameters parameters) : MooseApp(parameters)
{
  mece789App::registerAll(_factory, _action_factory, _syntax);
}

mece789App::~mece789App() {}

void
mece789App::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<mece789App>(f, af, syntax);
  Registry::registerObjectsTo(f, {"mece789App"});
  Registry::registerActionsTo(af, {"mece789App"});

  /* register custom execute flags, action syntax, etc. here */
}

void
mece789App::registerApps()
{
  registerApp(mece789App);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
mece789App__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  mece789App::registerAll(f, af, s);
}
extern "C" void
mece789App__registerApps()
{
  mece789App::registerApps();
}
