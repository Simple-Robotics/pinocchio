//
// Copyright (c) 2015-2021 CNRS INRIA
//

#ifndef __pinocchio_python_multibody_joint_joints_liegroup_hpp__
#define __pinocchio_python_multibody_joint_joints_liegroup_hpp__

#include <boost/python.hpp>

#include "pinocchio/bindings/python/fwd.hpp"
#include "pinocchio/multibody/liegroups.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    template<class JointModelDerived>
    struct JointModelLieGroupPythonVisitor
    : public boost::python::def_visitor<JointModelLieGroupPythonVisitor<JointModelDerived>>
    {
    public:
      typedef CartesianProductOperationVariantTpl<
        context::Scalar,
        context::Options,
        LieGroupCollectionDefaultTpl>
        LieGroupOperation;
      typedef typename LieGroupMap::template operation<JointModelDerived>::type LieGroupType;

      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl.def("liegroup", &liegroup).staticmethod("liegroup");
      }

      static LieGroupOperation liegroup()
      {
        return LieGroupOperation(LieGroupType());
      }
    };

  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_multibody_joint_joints_liegroup_hpp__
