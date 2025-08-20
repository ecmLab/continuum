#include "ecmApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
ecmApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  params.set<bool>("use_legacy_dirichlet_bc") = false;

  return params;
}

ecmApp::ecmApp(InputParameters parameters) : MooseApp(parameters)
{
  ecmApp::registerAll(_factory, _action_factory, _syntax);
}

ecmApp::~ecmApp() {}

void
ecmApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"ecmApp"});
  Registry::registerActionsTo(af, {"ecmApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
ecmApp::registerApps()
{
  registerApp(ecmApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
ecmApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ecmApp::registerAll(f, af, s);
}
extern "C" void
ecmApp__registerApps()
{
  ecmApp::registerApps();
}
