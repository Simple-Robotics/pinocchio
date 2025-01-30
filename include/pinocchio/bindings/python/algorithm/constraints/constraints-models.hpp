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

    template<class T>
    struct WrapperModelBase
    {
      typedef typename T::Base Base;
    }
    template<>
    struct WrapperModelBase<boost::blank>
    {
      typedef typename boost::blank Base;
    }

    // Add the inheritance
    template<class T>
    inline bp::class_<T> & expose_constraint_model_inheritance(bp::class_<T> & cl)
    {
      return ConstraintModelInheritanceHelper<T, typename WrapperModelBase<T>::Base>::expose_inheritance();
    }

    // generic expose_constraint_model : do nothing special
    template<class T>
    bp::class_<T> & expose_constraint_model(bp::class_<T> & cl)
    {
      return cl
      ;
    }

    // specialization for ConstraintModelRevolute
    // template<>
    // bp::class_<context::ConstraintModelRX> &
    // expose_constraint_model<context::ConstraintModelRX>(bp::class_<context::ConstraintModelRX> & cl)
    // {
    //   return cl
    //     .def(bp::init<>(
    //       bp::args("self"), "Init ConstraintModelRX with the X axis ([1, 0, 0]) as rotation axis."))
    //     .def(
    //       "getMotionAxis", &context::ConstraintModelRX::getMotionAxis,
    //       "Rotation axis of the ConstraintModelRX.");
    // };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_models_hpp__
