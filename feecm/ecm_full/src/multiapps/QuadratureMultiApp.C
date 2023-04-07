/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.cc to edit this template
 */

/* 
 * File:   QuadratureMultiApp.C
 * Author: srinathcs
 * 
 * Created on March 22, 2022, 2:34 PM
 */

//MOOSE includes
#include "QuadratureMultiApp.h"
#include "MooseMesh.h"
#include "FEProblem.h"

// libMesh includes
#include "libmesh/parallel_algebra.h"

registerMooseObject("MooseApp", QuadratureMultiApp);

InputParameters
QuadratureMultiApp::validParams()
{
  InputParameters params = TransientMultiApp::validParams();
  params += BlockRestrictable::validParams();
  params.addClassDescription(
      "Automatically generates Sub-App positions from centroids of elements in the master mesh.");
  params.suppressParameter<std::vector<Point>>("positions");
  params.suppressParameter<std::vector<FileName>>("positions_file");
  return params;
}

QuadratureMultiApp::QuadratureMultiApp(const InputParameters & parameters)
  : TransientMultiApp(parameters), BlockRestrictable(this)
{
}

void
QuadratureMultiApp::fillPositions()
{
  MooseMesh & master_mesh = _fe_problem.mesh();

  for (const auto & elem_ptr : master_mesh.getMesh().active_local_element_ptr_range())
  {
    if (hasBlocks(elem_ptr->subdomain_id()))
    {
        const FEFamily mapping_family = FEMap::map_fe_type(*elem_ptr);
        FEType fe_type(elem_ptr->default_order(),mapping_family);
        std::unique_ptr<FEBase> fe(FEBase::build(elem_ptr->dim(), fe_type));
        const std::vector<Point> & q_points = fe->get_xyz();
        const int extraorder = 0;
        std::unique_ptr<QBase> qrule (fe_type.default_quadrature_rule (2, extraorder));
        fe->attach_quadrature_rule (qrule.get());
        fe->reinit(elem_ptr);
        for (auto i=0; i<q_points.size(); ++i){ 
            _positions.push_back(q_points[i]);
        }
    }
      
  }
  // Use the comm from the problem this MultiApp is part of
  libMesh::ParallelObject::comm().allgather(_positions);

  if (_positions.empty())
    mooseError("No positions found for QuadratureMultiapp ", _name);

  // An attempt to try to make this parallel stable
  std::sort(_positions.begin(), _positions.end());
}
