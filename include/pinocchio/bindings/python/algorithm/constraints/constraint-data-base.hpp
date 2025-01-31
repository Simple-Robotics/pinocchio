//
// Copyright (c) 2015-2025 CNRS INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_data_base_hpp__
#define __pinocchio_python_algorithm_constraints_data_base_hpp__

#include <boost/python.hpp>
#include <eigenpy/exception.hpp>
#include <eigenpy/eigen-to-python.hpp>

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/constraint-data-base.hpp"
#include "pinocchio/bindings/python/fwd.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    template<class ConstraintDataDerived>
    struct ConstraintDataBasePythonVisitor
    : public bp::def_visitor<ConstraintDataBasePythonVisitor<ConstraintDataDerived>>
    {
      typedef ConstraintDataDerived Self;
    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl
          .def("classname", &Self::classname)
          .staticmethod("classname")
          .def("shortname", &Self::shortname)
#ifndef PINOCCHIO_PYTHON_SKIP_COMPARISON_OPERATIONS
          .def(bp::self == bp::self)
          .def(bp::self != bp::self)
#endif
          ;
      }
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_data_base_hpp__
