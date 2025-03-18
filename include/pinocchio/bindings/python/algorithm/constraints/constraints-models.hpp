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
          bp::return_value_policy<bp::copy_const_reference>());
      return cl;
    }

    template<>
    bp::class_<context::JointLimitConstraintModel> &
    expose_constraint_model(bp::class_<context::JointLimitConstraintModel> & cl)
    {
      typedef typename context::JointLimitConstraintModel::JointIndexVector JointIndexVector;
      typedef typename context::JointLimitConstraintModel Self;
      cl.def(bp::init<const context::Model &, const JointIndexVector &>(
               (bp::arg("self"), bp::arg("model"), bp::arg("activable_joints")),
               "Contructor from given joint index vector "
               "implied in the constraint."))
        .def(
          "getActiveSetIndexes", &Self::getActiveSetIndexes,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Indexes of the active constraints set.")
        .def(
          "getActivableLowerBoundConstraints", &Self::getActivableLowerBoundConstraints,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Returns the vector of configuration vector index for activable lower bounds.")
        .def(
          "getActiveLowerBoundConstraints", &Self::getActiveLowerBoundConstraints,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Returns the vector of configuration vector index for active lower bounds.")
        .def(
          "getActivableUpperBoundConstraints", &Self::getActivableUpperBoundConstraints,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Returns the vector of configuration vector index for activable upper bounds.")
        .def(
          "getActiveUpperBoundConstraints", &Self::getActiveUpperBoundConstraints,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Returns the vector of configuration vector index for active upper bounds.")
        .def(
          "getActivableLowerBoundConstraintsTangent",
          &Self::getActivableLowerBoundConstraintsTangent,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Returns the vector of tangent configuration vector index for activable lower bounds.")
        .def(
          "getActiveLowerBoundConstraintsTangent", &Self::getActiveLowerBoundConstraintsTangent,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Returns the vector of tangent configuration vector index for active lower bounds.")
        .def(
          "getActivableUpperBoundConstraintsTangent",
          &Self::getActivableUpperBoundConstraintsTangent,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Returns the vector of tangent configuration vector index for activable upper bounds.")
        .def(
          "getActiveUpperBoundConstraintsTangent", &Self::getActiveUpperBoundConstraintsTangent,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Returns the vector of tangent configuration vector index for active upper bounds.");
      // resize
      return cl;
    }
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_models_hpp__
