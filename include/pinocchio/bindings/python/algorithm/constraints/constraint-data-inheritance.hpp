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

    // Default inheritance Visitor Template
    template<class ConstraintDataDerived, class ConstraintDataBase>
    struct ConstraintDataInheritancePythonVisitor
    : public bp::def_visitor<ConstraintDataInheritancePythonVisitor<ConstraintDataDerived, ConstraintDataBase>>
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
    template<class ConstraintDataDerived>
    struct ConstraintDataInheritancePythonVisitor<ConstraintDataDerived, FrameConstraintDataBase<ConstraintDataDerived>>
    {
    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl
          ;
      }
    };

    template<class ConstraintDataDerived>
    struct ConstraintDataInheritancePythonVisitor<ConstraintDataDerived, PointConstraintDataBase<ConstraintDataDerived>>
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

#endif // ifndef __pinocchio_python_algorithm_constraints_data_inheritance_hpp__
