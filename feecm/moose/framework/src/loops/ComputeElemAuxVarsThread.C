//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeElemAuxVarsThread.h"

// MOOSE includes
#include "AuxiliarySystem.h"
#include "AuxKernel.h"
#include "SwapBackSentinel.h"
#include "FEProblem.h"
#include "MaterialBase.h"

#include "libmesh/threads.h"

template <typename AuxKernelType>
Threads::spin_mutex ComputeElemAuxVarsThread<AuxKernelType>::writable_variable_mutex;

template <typename AuxKernelType>
ComputeElemAuxVarsThread<AuxKernelType>::ComputeElemAuxVarsThread(
    FEProblemBase & problem,
    const MooseObjectWarehouse<AuxKernelType> & storage,
    const std::vector<std::vector<MooseVariableFEBase *>> & vars,
    bool need_materials)
  : ThreadedElementLoop<ConstElemRange>(problem),
    _aux_sys(problem.getAuxiliarySystem()),
    _aux_kernels(storage),
    _aux_vars(vars),
    _need_materials(need_materials)
{
}

// Splitting Constructor
template <typename AuxKernelType>
ComputeElemAuxVarsThread<AuxKernelType>::ComputeElemAuxVarsThread(ComputeElemAuxVarsThread & x,
                                                                  Threads::split /*split*/)
  : ThreadedElementLoop<ConstElemRange>(x._fe_problem),
    _aux_sys(x._aux_sys),
    _aux_kernels(x._aux_kernels),
    _aux_vars(x._aux_vars),
    _need_materials(x._need_materials)
{
}

template <typename AuxKernelType>
ComputeElemAuxVarsThread<AuxKernelType>::~ComputeElemAuxVarsThread()
{
}

template <typename AuxKernelType>
void
ComputeElemAuxVarsThread<AuxKernelType>::subdomainChanged()
{
  _fe_problem.subdomainSetup(_subdomain, _tid);

  // prepare variables
  for (auto * var : _aux_vars[_tid])
    var->prepareAux();

  std::set<MooseVariableFEBase *> needed_moose_vars;
  std::set<unsigned int> needed_mat_props;
  std::set<TagID> needed_fe_var_matrix_tags;
  std::set<TagID> needed_fe_var_vector_tags;

  _fe_problem.getMaterialWarehouse().updateBlockFEVariableCoupledVectorTagDependency(
      _subdomain, needed_fe_var_vector_tags, _tid);
  if (_aux_kernels.hasActiveBlockObjects(_subdomain, _tid))
  {
    const std::vector<std::shared_ptr<AuxKernelType>> & kernels =
        _aux_kernels.getActiveBlockObjects(_subdomain, _tid);
    for (const auto & aux : kernels)
    {
      aux->subdomainSetup();
      const std::set<MooseVariableFEBase *> & mv_deps = aux->getMooseVariableDependencies();
      const std::set<unsigned int> & mp_deps = aux->getMatPropDependencies();
      needed_moose_vars.insert(mv_deps.begin(), mv_deps.end());
      needed_mat_props.insert(mp_deps.begin(), mp_deps.end());

      auto & fe_var_coup_vtags = aux->getFEVariableCoupleableVectorTags();
      needed_fe_var_vector_tags.insert(fe_var_coup_vtags.begin(), fe_var_coup_vtags.end());

      auto & fe_var_coup_mtags = aux->getFEVariableCoupleableMatrixTags();
      needed_fe_var_matrix_tags.insert(fe_var_coup_mtags.begin(), fe_var_coup_mtags.end());
    }
  }

  _fe_problem.setActiveElementalMooseVariables(needed_moose_vars, _tid);
  _fe_problem.setActiveMaterialProperties(needed_mat_props, _tid);
  _fe_problem.prepareMaterials(_subdomain, _tid);
  _fe_problem.setActiveFEVariableCoupleableMatrixTags(needed_fe_var_matrix_tags, _tid);
  _fe_problem.setActiveFEVariableCoupleableVectorTags(needed_fe_var_vector_tags, _tid);
}

template <typename AuxKernelType>
void
ComputeElemAuxVarsThread<AuxKernelType>::onElement(const Elem * elem)
{
  if (_aux_kernels.hasActiveBlockObjects(_subdomain, _tid))
  {
    const std::vector<std::shared_ptr<AuxKernelType>> & kernels =
        _aux_kernels.getActiveBlockObjects(_subdomain, _tid);
    _fe_problem.prepare(elem, _tid);
    _fe_problem.reinitElem(elem, _tid);

    // Set up the sentinel so that, even if reinitMaterials() throws, we
    // still remember to swap back.
    SwapBackSentinel sentinel(_fe_problem, &FEProblem::swapBackMaterials, _tid, _need_materials);

    if (_need_materials)
      _fe_problem.reinitMaterials(elem->subdomain_id(), _tid);

    for (const auto & aux : kernels)
    {
      aux->compute();

      // update the aux solution vector if writable coupled variables are used
      if (aux->hasWritableCoupledVariables())
      {
        Threads::spin_mutex::scoped_lock lock(writable_variable_mutex);
        for (auto * var : aux->getWritableCoupledVariables())
        {
          // insert into the global solution vector
          var->insert(_aux_sys.solution());
          var->prepareAux();
        }

        // make solution values available for dependent AuxKernels
        aux->variable().insert(_aux_sys.solution());
        aux->variable().prepareAux();
        _fe_problem.reinitElem(elem, _tid);
      }
    }

    // update the solution vector
    {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      for (auto * var : _aux_vars[_tid])
        var->insert(_aux_sys.solution());
    }
  }
}

template <typename AuxKernelType>
void
ComputeElemAuxVarsThread<AuxKernelType>::post()
{
  _fe_problem.clearActiveElementalMooseVariables(_tid);
  _fe_problem.clearActiveMaterialProperties(_tid);

  _fe_problem.clearActiveFEVariableCoupleableVectorTags(_tid);
  _fe_problem.clearActiveFEVariableCoupleableMatrixTags(_tid);
}

template <typename AuxKernelType>
void
ComputeElemAuxVarsThread<AuxKernelType>::join(const ComputeElemAuxVarsThread & /*y*/)
{
}

template class ComputeElemAuxVarsThread<AuxKernel>;
template class ComputeElemAuxVarsThread<VectorAuxKernel>;
template class ComputeElemAuxVarsThread<ArrayAuxKernel>;
