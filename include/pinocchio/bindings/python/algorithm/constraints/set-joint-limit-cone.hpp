//
// Copyright (c) 2022 INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_set_joint_limit_cone_hpp__
#define __pinocchio_python_algorithm_constraints_set_joint_limit_cone_hpp__

#include "pinocchio/algorithm/constraints/joint-limit-constraint-cone.hpp"

#include "pinocchio/bindings/python/algorithm/constraints/set-base.hpp"
#include "pinocchio/bindings/python/utils/cast.hpp"
#include "pinocchio/bindings/python/utils/copyable.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    template<typename JointLimitConstraintCone>
    struct JointLimitConstraintConePythonVisitor
    : public boost::python::def_visitor<
        JointLimitConstraintConePythonVisitor<JointLimitConstraintCone>>
    {

      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl.def(bp::init<Eigen::DenseIndex, Eigen::DenseIndex>(
          bp::args("self", "negative_orthant_size", "positive_orthant_size"),
          "Default constructor given a positive and a negative size."));
        // resize
        // conservativeresize
      }

      static void expose()
      {
        bp::class_<JointLimitConstraintCone>(
          "JointLimitConstraintCone", "Concatenation of positive and negative orthant cone",
          bp::no_init)
          .def(SetPythonVisitor<JointLimitConstraintCone, context::VectorXs>())
          .def(ConeSetPythonVisitor<JointLimitConstraintCone>())
          .def(JointLimitConstraintConePythonVisitor())
          // .def(CastVisitor<JointLimitConstraintCone>())
          // .def(ExposeConstructorByCastVisitor<JointLimitConstraintCone,::pinocchio::JointLimitConstraintCone>())
          .def(CopyableVisitor<JointLimitConstraintCone>());
      }
    };

  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_set_joint_limit_cone_hpp__
