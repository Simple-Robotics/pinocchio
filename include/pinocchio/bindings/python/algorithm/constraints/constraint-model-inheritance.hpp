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

    // Default inheritance Visitor Template
    template<class ConstraintModelDerived, class ConstraintModelBase>
    struct ConstraintModelInheritancePythonVisitor
    : public bp::def_visitor<ConstraintModelInheritancePythonVisitor<ConstraintModelDerived, ConstraintModelBase>>
    {
    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl
          ;
      }
    };

    // Specialize
    template<class ConstraintModelDerived>
    struct ConstraintModelInheritancePythonVisitor<ConstraintModelDerived, FrameConstraintModelBase<ConstraintModelDerived>>
    {
    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl
          ;
      }
    };

    template<class ConstraintModelDerived>
    struct ConstraintModelInheritancePythonVisitor<ConstraintModelDerived, PointConstraintModelBase<ConstraintModelDerived>>
    {
    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl
          ;
      }
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_model_inheritance_hpp__
