//
// Copyright (c) 2015-2025 CNRS INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_model_inheritance_hpp__
#define __pinocchio_python_algorithm_constraints_model_inheritance_hpp__

#include <boost/python.hpp>
#include <eigenpy/exception.hpp>
#include <eigenpy/eigen-to-python.hpp>

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/frame-constraint-model-base.hpp"
#include "pinocchio/algorithm/constraints/point-constraint-model-base.hpp"
#include "pinocchio/bindings/python/fwd.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    // Default inheritance template
    template<class ConstraintModelDerived, class ConstraintModelBase>
    struct ConstraintModelInheritanceHelper
    {
      static bp::class_<ConstraintModelDerived> & expose_inheritance(bp::class_<ConstraintModelDerived> & cl)
      {
        return cl
        ;
      }
    };

    // Specialize
    template<class ConstraintModelDerived>
    struct ConstraintModelInheritanceHelper<ConstraintModelDerived, FrameConstraintModelBase<ConstraintModelDerived>>
    {
      static bp::class_<ConstraintModelDerived> & expose_inheritance(bp::class_<ConstraintModelDerived> & cl)
      {
        return cl
        ;
      }
    };

    template<class ConstraintModelDerived>
    struct ConstraintModelInheritanceHelper<ConstraintModelDerived, PointConstraintModelBase<ConstraintModelDerived>>
    {
      static bp::class_<ConstraintModelDerived> & expose_inheritance(bp::class_<ConstraintModelDerived> & cl)
      {
        return cl
        ;
      }
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_model_inheritance_hpp__
