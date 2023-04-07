#include "ecmElectrochemApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
ecmElectrochemApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy material output, i.e., output properties on INITIAL as well as TIMESTEP_END
  params.set<bool>("use_legacy_material_output") = false;

  return params;
}

ecmElectrochemApp::ecmElectrochemApp(InputParameters parameters) : MooseApp(parameters)
{
  ecmElectrochemApp::registerAll(_factory, _action_factory, _syntax);
}

ecmElectrochemApp::~ecmElectrochemApp() {}

void
ecmElectrochemApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"ecmElectrochemApp"});
  Registry::registerActionsTo(af, {"ecmElectrochemApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
ecmElectrochemApp::registerApps()
{
  registerApp(ecmElectrochemApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
ecmElectrochemApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ecmElectrochemApp::registerAll(f, af, s);
}
extern "C" void
ecmElectrochemApp__registerApps()
{
  ecmElectrochemApp::registerApps();
}
