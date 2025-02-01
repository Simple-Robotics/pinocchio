//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_datas_hpp__
#define __pinocchio_python_algorithm_constraints_datas_hpp__

#include "pinocchio/algorithm/constraints/constraint-data-generic.hpp"
#include "pinocchio/algorithm/constraints/constraint-collection-default.hpp"
#include "pinocchio/bindings/python/fwd.hpp"
#include "pinocchio/bindings/python/algorithm/constraints/constraint-data-inheritance.hpp"
#include "pinocchio/bindings/python/utils/printable.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    // generic expose_constraint_data : do nothing special
    template<class T>
    inline bp::class_<T> & expose_constraint_data(bp::class_<T> & cl)
    {
      return cl;
    }

    // specialization for ConstraintDatas
    template<>
    bp::class_<context::BilateralPointConstraintData> &
    expose_constraint_data(bp::class_<context::BilateralPointConstraintData> & cl)
    {
      return cl;
    }

    template<>
    bp::class_<context::FrictionalPointConstraintData> &
    expose_constraint_data(bp::class_<context::FrictionalPointConstraintData> & cl)
    {
      return cl;
    }

    template<>
    bp::class_<context::WeldConstraintData> &
    expose_constraint_data(bp::class_<context::WeldConstraintData> & cl)
    {
      return cl;
    }

    template<>
    bp::class_<context::FrictionalJointConstraintData> &
    expose_constraint_data(bp::class_<context::FrictionalJointConstraintData> & cl)
    {
      return cl.def(
        bp::init<const typename context::FrictionalJointConstraintData::ConstraintModel &>(
          bp::args("self", "constraint_model")));
    }

    template<>
    bp::class_<context::JointLimitConstraintData> &
    expose_constraint_data(bp::class_<context::JointLimitConstraintData> & cl)
    {
      return cl
        .def(bp::init<const typename context::JointLimitConstraintData::ConstraintModel &>(
          bp::args("self", "constraint_model")))
        .PINOCCHIO_ADD_PROPERTY(
          context::JointLimitConstraintData, constraint_residual, "Constraint residual.");
    }
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_datas_hpp__
