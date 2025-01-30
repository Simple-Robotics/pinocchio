//
// Copyright (c) 2015-2025 CNRS INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_model_base_hpp__
#define __pinocchio_python_algorithm_constraints_model_base_hpp__

#include <boost/python.hpp>
#include <eigenpy/exception.hpp>
#include <eigenpy/eigen-to-python.hpp>

#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-base.hpp"
#include "pinocchio/bindings/python/fwd.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    template<class ConstraintModelDerived>
    struct ConstraintModelBasePythonVisitor
    : public bp::def_visitor<ConstraintModelBasePythonVisitor<ConstraintModelDerived>>
    {
    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl
          // .def(bp::init<>(bp::arg("self")))
          // // All are add_properties cause ReadOnly
          // .add_property("id", &get_id)
          // .add_property(
          //   "hasConfigurationLimit", &ConstraintModelDerived::hasConfigurationLimit,
          //   "Return vector of boolean if joint has configuration limits.")
          // .def("classname", &ConstraintModelDerived::classname)
          // .staticmethod("classname")
          // .def("calc", &calc0, bp::args("self", "cdata", "q"))
#ifndef PINOCCHIO_PYTHON_SKIP_COMPARISON_OPERATIONS
          .def(bp::self == bp::self)
          .def(bp::self != bp::self)
#endif
          ;
      // }
      // static JointIndex get_id(const ConstraintModelDerived & self)
      // {
      //   return self.id();
      // }
      // static void
      // calc0(const ConstraintModelDerived & self, ConstraintDataDerived & cdata, const context::VectorXs & q)
      // {
      //   self.calc(cdata, q);
      // }
      // static void calc1(
      //   const ConstraintModelDerived & self,
      //   ConstraintDataDerived & cdata,
      //   const context::VectorXs & q,
      //   const context::VectorXs & v)
      // {
      //   self.calc(cdata, q, v);
      }
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_model_base_hpp__
