//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_models_hpp__
#define __pinocchio_python_algorithm_constraints_models_hpp__

#include "pinocchio/algorithm/constraints/constraint-model-generic.hpp"
#include "pinocchio/algorithm/constraints/constraint-collection-default.hpp"
#include "pinocchio/bindings/python/fwd.hpp"
#include "pinocchio/bindings/python/utils/printable.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    // generic expose_constraint_model : do nothing special
    template<class T>
    bp::class_<T> & expose_constraint_model(bp::class_<T> & cl)
    {
      return cl;
    }

    // specialization for ConstraintModels
    template<>
    bp::class_<context::BilateralPointConstraintModel> &
    expose_constraint_model(bp::class_<context::BilateralPointConstraintModel> & cl)
    {
      return cl;
    }

    template<>
    bp::class_<context::FrictionalPointConstraintModel> &
    expose_constraint_model(bp::class_<context::FrictionalPointConstraintModel> & cl)
    {
      return cl;
    }

    template<>
    bp::class_<context::WeldConstraintModel> &
    expose_constraint_model(bp::class_<context::WeldConstraintModel> & cl)
    {
      return cl;
    }

    template<>
    bp::class_<context::FrictionalJointConstraintModel> &
    expose_constraint_model(bp::class_<context::FrictionalJointConstraintModel> & cl)
    {

      typedef typename context::FrictionalJointConstraintModel::JointIndexVector JointIndexVector;

      cl.def(bp::init<const context::Model &, const JointIndexVector &>(
               (bp::arg("self"), bp::arg("model"), bp::arg("active_joints")),
               "Contructor from given joint index vector "
               "implied in the constraint."))
        .def(
          "getActiveDofs", &context::FrictionalJointConstraintModel::getActiveDofs,
          bp::return_value_policy<bp::copy_const_reference>())
        // .def("getRowActiveIndexes", ...)
        // .def("getRowSparsityPattern", ...)
        ;
      return cl;
    }

    template<>
    bp::class_<context::JointLimitConstraintModel> &
    expose_constraint_model(bp::class_<context::JointLimitConstraintModel> & cl)
    {
      typedef typename context::FrictionalJointConstraintModel::JointIndexVector JointIndexVector;
      cl.def(bp::init<const context::Model &, const JointIndexVector &>(
               (bp::arg("self"), bp::arg("model"), bp::arg("active_joints")),
               "Contructor from given joint index vector "
               "implied in the constraint."))
        // Init with lb and ub
        .def(
          "getActiveLowerBoundConstraints",
          &context::JointLimitConstraintModel::getActiveLowerBoundConstraints,
          bp::return_value_policy<bp::copy_const_reference>(), "Active lower bound constraints.")
        .def(
          "getActiveLowerBoundConstraintsTangent",
          &context::JointLimitConstraintModel::getActiveLowerBoundConstraintsTangent,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Active lower bound constraints in tangent.")
        .def(
          "getActiveUpperBoundConstraints",
          &context::JointLimitConstraintModel::getActiveUpperBoundConstraints,
          bp::return_value_policy<bp::copy_const_reference>(), "Active upper bound constraints.")
        .def(
          "getActiveUpperBoundConstraintsTangent",
          &context::JointLimitConstraintModel::getActiveUpperBoundConstraintsTangent,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Active upper bound constraints in tangent.")
        // .def("getRowActiveIndexes", ...)
        // .def("getRowSparsityPattern", ...)
        ;
      return cl;
    }
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_models_hpp__
