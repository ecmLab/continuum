
#include "CZMDiffusiveVariableKernelBase.h"

InputParameters
CZMDiffusiveVariableKernelBase::validParams()
{
    InputParameters params = InterfaceKernel::validParams();
      params.addRequiredParam<unsigned int>("component",
                                        "the component of the "
                                        "displacement vector this kernel is working on:"
                                        " component == 0, ==> X"
                                        " component == 1, ==> Y"
                                        " component == 2, ==> Z");
    params.set<bool>("_use_undisplaced_reference_points") = true;
//    params.addRequiredCoupledVar("var", "Name of variable to couple czm to");
    params.addCoupledVar("displacements", "the string containing displacement variables");
    params.addParam<std::string>("base_name", "Material property base name");
    params.set<std::string>("flux_global_name") = "flux_global";
    params.addParam<bool>("include_gap", true, "include displacement"
        "contribution to the off diagonal jacobian, setting to false will only "
        "include diagonal jacobian components");
    return params;
}

CZMDiffusiveVariableKernelBase::CZMDiffusiveVariableKernelBase(const InputParameters& parameters)
    : JvarMapKernelInterface<InterfaceKernel>(parameters),
    _base_name(isParamValid("base_name") && !getParam<std::string>("base_name").empty()
                   ? getParam<std::string>("base_name") + "_"
                   : ""),
    _include_gap(getParam<bool>("include_gap")),
    _component(_include_gap ? getParam<unsigned int>("component") : 0),
    _ndisp(_include_gap ? coupledComponents("displacements") : 0),
    _disp_var(_ndisp),
    _disp_neighbor_var(_ndisp),
    _disp_vars(_ndisp),
    _flux_global(getMaterialPropertyByName<Real>(
        _base_name + getParam<std::string>("flux_global_name"))),
    _dflux_dvariablejump_global(getMaterialPropertyByName<Real>(
        _base_name +"dflux_dvariablejump_global")),
    _dflux_djump_global(
        getMaterialPropertyByName<RealVectorValue>(_base_name + "dflux_djump_global"))
//    _var_id(coupled("variable"))
{
  if (_include_gap)
  {
    if (_ndisp != _mesh.dimension())
      paramError("displacements", "Number of displacements must match problem dimension.");

    if (_ndisp > 3 || _ndisp < 1)
      mooseError("the CZM material requires 1, 2 or 3 displacement variables");

    for (unsigned int i = 0; i < _ndisp; ++i)
    {
      _disp_var[i] = coupled("displacements", i);
      _disp_neighbor_var[i] = coupled("displacements", i);
      _disp_vars[i] = getVar("displacements", i);
    }
  }
//    _var = getVar("var", 0);
}

Real
CZMDiffusiveVariableKernelBase::computeQpResidual(Moose::DGResidualType type)
{
  // For diffusive kernels, the flux term is a scalar measure
  Real r = _flux_global[_qp];
  switch (type)
  {
    // [test_secondary-test_primary]*T where T represents the traction.
    case Moose::Element:
      r *= -_test[_i][_qp];
      break;
    case Moose::Neighbor:
      r *= _test_neighbor[_i][_qp];
      break;
  }
  return r;
}

Real
CZMDiffusiveVariableKernelBase::computeQpJacobian(Moose::DGJacobianType type)
{
  // retrieve the diagonal Jacobian coefficient dependning value of local variable
  // jump
  return computeDResidualDVariable(type);
}

Real
CZMDiffusiveVariableKernelBase::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
  // bail out if jvar is not coupled
  if (getJvarMap()[jvar] < 0)
    return 0.0;
  if (_include_gap)
  {
    // Jacobian of the residual w.r.t to the coupled displacement
    // component
    for (unsigned int off_diag_component = 0; off_diag_component < _ndisp; ++off_diag_component)
    {
      if (jvar == _disp_var[off_diag_component])
        return computeDResidualDDisplacement(off_diag_component, type);
    }
  }
  // this is the place where one should implement derivatives of the residual w.r.t. other variables
  return 0.0;
}
