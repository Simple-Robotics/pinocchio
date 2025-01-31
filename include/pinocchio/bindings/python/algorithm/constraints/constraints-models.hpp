//
// Copyright (c) 2015-2025 CNRS INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_models_hpp__
#define __pinocchio_python_algorithm_constraints_models_hpp__

#include "pinocchio/algorithm/constraints/constraint-model-generic.hpp"
#include "pinocchio/algorithm/constraints/constraint-collection-default.hpp"
#include "pinocchio/bindings/python/algorithm/constraints/constraint-model-inheritance.hpp"
#include "pinocchio/bindings/python/utils/printable.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    // Add the inheritance
    template<class T>
    inline bp::class_<T> & expose_constraint_model_inheritance(bp::class_<T> & cl)
    {
      return cl
      .def(ConstraintModelInheritancePythonVisitor<T, typename T::Base>());
      ;
    }

    // generic expose_constraint_model : do nothing special
    template<class T>
    bp::class_<T> & expose_constraint_model(bp::class_<T> & cl)
    {
      return cl
      ;
    }

    // specialization for ConstraintModels
    template<>
    bp::class_<context::BilateralPointConstraintModel> & expose_constraint_model(
      bp::class_<context::BilateralPointConstraintModel> & cl)
    {
      return cl
      ;
    }

    template<>
    bp::class_<context::FrictionalPointConstraintModel> & expose_constraint_model(
      bp::class_<context::FrictionalPointConstraintModel> & cl)
    {
      return cl
      ;
    }

    template<>
    bp::class_<context::FrictionalJointConstraintModel> & expose_constraint_model(
      bp::class_<context::FrictionalJointConstraintModel> & cl)
    {
      return cl
      ;
    }

    template<>
    bp::class_<context::JointLimitConstraintModel> & expose_constraint_model(
      bp::class_<context::JointLimitConstraintModel> & cl)
    {
      return cl
      ;
    }

    template<>
    bp::class_<context::WeldConstraintModel> & expose_constraint_model(
      bp::class_<context::WeldConstraintModel> & cl
    )
    {
      return cl
      ;
    }
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_models_hpp__
