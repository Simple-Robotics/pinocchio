//
// Copyright (c) 2015-2025 CNRS INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_data_inheritance_hpp__
#define __pinocchio_python_algorithm_constraints_data_inheritance_hpp__

#include <boost/python.hpp>
#include <eigenpy/exception.hpp>
#include <eigenpy/eigen-to-python.hpp>

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/frame-constraint-data-base.hpp"
#include "pinocchio/algorithm/constraints/point-constraint-data-base.hpp"
#include "pinocchio/bindings/python/fwd.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    // Default inheritance template
    template<class ConstraintDataDerived, class ConstraintDataBase>
    struct ConstraintDataInheritanceHelper
    {
      static bp::class_<ConstraintDataDerived> & expose_inheritance(bp::class_<ConstraintDataDerived> & cl)
      {
        return cl
        ;
      }
    };

    // Specialize
    template<class ConstraintDataDerived>
    struct ConstraintDataInheritanceHelper<ConstraintDataDerived, FrameConstraintDataBase<ConstraintDataDerived>>
    {
      static bp::class_<ConstraintDataDerived> & expose_inheritance(bp::class_<ConstraintDataDerived> & cl)
      {
        return cl
        ;
      }
    };

    template<class ConstraintDataDerived>
    struct ConstraintDataInheritanceHelper<ConstraintDataDerived, PointConstraintDataBase<ConstraintDataDerived>>
    {
      static bp::class_<ConstraintDataDerived> & expose_inheritance(bp::class_<ConstraintDataDerived> & cl)
      {
        return cl
        ;
      }
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_data_inheritance_hpp__
