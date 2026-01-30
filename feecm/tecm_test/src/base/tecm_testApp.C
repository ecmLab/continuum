#include "tecm_testApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
tecm_testApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

tecm_testApp::tecm_testApp(InputParameters parameters) : MooseApp(parameters)
{
  tecm_testApp::registerAll(_factory, _action_factory, _syntax);
}

tecm_testApp::~tecm_testApp() {}

void
tecm_testApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<tecm_testApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"tecm_testApp"});
  Registry::registerActionsTo(af, {"tecm_testApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
tecm_testApp::registerApps()
{
  registerApp(tecm_testApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
tecm_testApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  tecm_testApp::registerAll(f, af, s);
}
extern "C" void
tecm_testApp__registerApps()
{
  tecm_testApp::registerApps();
}
