#include "edl_solidApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
edl_solidApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

edl_solidApp::edl_solidApp(InputParameters parameters) : MooseApp(parameters)
{
  edl_solidApp::registerAll(_factory, _action_factory, _syntax);
}

edl_solidApp::~edl_solidApp() {}

void 
edl_solidApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<edl_solidApp>(f, af, s);
  Registry::registerObjectsTo(f, {"edl_solidApp"});
  Registry::registerActionsTo(af, {"edl_solidApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
edl_solidApp::registerApps()
{
  registerApp(edl_solidApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
edl_solidApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  edl_solidApp::registerAll(f, af, s);
}
extern "C" void
edl_solidApp__registerApps()
{
  edl_solidApp::registerApps();
}
