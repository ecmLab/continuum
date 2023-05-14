#include "ecm_fullApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
ecm_fullApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  params.set<bool>("use_legacy_dirichlet_bc") = false;

  return params;
}

ecm_fullApp::ecm_fullApp(InputParameters parameters) : MooseApp(parameters)
{
  ecm_fullApp::registerAll(_factory, _action_factory, _syntax);
}

ecm_fullApp::~ecm_fullApp() {}

void
ecm_fullApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"ecm_fullApp"});
  Registry::registerActionsTo(af, {"ecm_fullApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
ecm_fullApp::registerApps()
{
  registerApp(ecm_fullApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
ecm_fullApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ecm_fullApp::registerAll(f, af, s);
}
extern "C" void
ecm_fullApp__registerApps()
{
  ecm_fullApp::registerApps();
}
