#include "ecBetaApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
ecBetaApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  return params;
}

ecBetaApp::ecBetaApp(InputParameters parameters) : MooseApp(parameters)
{
  ecBetaApp::registerAll(_factory, _action_factory, _syntax);
}

ecBetaApp::~ecBetaApp() {}

void
ecBetaApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"ecBetaApp"});
  Registry::registerActionsTo(af, {"ecBetaApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
ecBetaApp::registerApps()
{
  registerApp(ecBetaApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
ecBetaApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ecBetaApp::registerAll(f, af, s);
}
extern "C" void
ecBetaApp__registerApps()
{
  ecBetaApp::registerApps();
}
