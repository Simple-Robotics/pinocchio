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
      typedef typename context::JointLimitConstraintModel::ConstraintData ConstraintData;
      typedef typename context::JointLimitConstraintModel Self;
      cl.def(bp::init<const context::Model &, const JointIndexVector &>(
               (bp::arg("self"), bp::arg("model"), bp::arg("activable_joints")),
               "Contructor from given joint index vector "
               "implied in the constraint."))
        .def("getNqReduce", &Self::getNqReduce, "Sum of nq of activable joints.")
        .def("getNvMaxAtom", &Self::getNvMaxAtom, "Max nv of atomic joints in activable joints.")
        .def("lowerSize", &Self::lowerSize, "Part of size() that are lower bound limits.")
        .def(
          "lowerActiveSize", &Self::lowerActiveSize,
          "Part of activeSize() that are lower bound limits.")
        .def("upperSize", &Self::upperSize, "Part of size() that are upper bound limits.")
        .def(
          "upperActiveSize", &Self::upperActiveSize,
          "Part of activeSize() that are upper bound limits.")
        .def(
          "getBoundPositionLimit", &Self::getBoundPositionLimit,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Position limit of the dof of the constraints.")
        .def(
          "getBoundPositionMargin", &Self::getBoundPositionMargin,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Position margin of the dof of the constraints.")
        .def(
          "getActivableJoints", &Self::getActivableJoints,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Joints for which there is at least one position limit.")
        .def(
          "getActivableIdxQs", &Self::getActivableIdxQs,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Q index in configuration of each limit constraint.")
        .def(
          "getActivableIdxQsReduce", &Self::getActivableIdxQsReduce,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Q index in thre reduce configuration (about activable joints) of each activable limit "
          "constraint.")
        .def(
          "getActiveIdxQsReduce", &Self::getActiveIdxQsReduce,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Q index in thre reduce configuration (about activable joints) of each active limit "
          "constraint.")
        .def(
          "getActivableNvs", &Self::getActivableNvs,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Nv of the atomic joint for which each activable position limit contribute to.")
        .def(
          "getActiveNvs", &Self::getActiveNvs, bp::return_value_policy<bp::copy_const_reference>(),
          "Nv of the atomic joint for which each active position limit contribute to.")
        .def(
          "getActivableIdxVs", &Self::getActivableIdxVs,
          bp::return_value_policy<bp::copy_const_reference>(),
          "V index of the atomic joint for which each activable position limit contribute to.")
        .def(
          "getActiveIdxVs", &Self::getActiveIdxVs,
          bp::return_value_policy<bp::copy_const_reference>(),
          "V index of the atomic joint for which each active position limit contribute to.")
        .def(
          "getActiveSetIndexes", &Self::getActiveSetIndexes,
          bp::return_value_policy<bp::copy_const_reference>(),
          "Indexes of the active constraints set.");
      return cl;
    }
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_models_hpp__
