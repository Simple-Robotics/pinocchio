//
// Copyright (c) 2022 INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_set_trivial_cones_hpp__
#define __pinocchio_python_algorithm_constraints_set_trivial_cones_hpp__

#include "pinocchio/algorithm/constraints/null-set.hpp"
#include "pinocchio/algorithm/constraints/unbounded-set.hpp"
#include "pinocchio/algorithm/constraints/orthant-cone.hpp"
#include "pinocchio/algorithm/constraints/joint-limit-constraint-cone.hpp"

#include "pinocchio/bindings/python/algorithm/constraints/set-base.hpp"
#include "pinocchio/bindings/python/utils/cast.hpp"
#include "pinocchio/bindings/python/utils/copyable.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    template<typename TrivialCone>
    struct TrivialConePythonVisitor
    : public boost::python::def_visitor<TrivialConePythonVisitor<TrivialCone>>
    {

      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl.def(bp::init<Eigen::DenseIndex>(bp::args("self", "size"), "Default constructor."));
        // resize
        // conservativeResize
      }

      static void expose(const std::string & class_name, const std::string & doc_string = "")
      {
        bp::class_<TrivialCone>(class_name.c_str(), doc_string.c_str(), bp::no_init)
          .def(SetPythonVisitor<TrivialCone, context::VectorXs>())
          .def(ConeSetPythonVisitor<TrivialCone>())
          .def(TrivialConePythonVisitor())
          // .def(CastVisitor<TrivialCone>())
          // .def(ExposeConstructorByCastVisitor<TrivialCone,::pinocchio::TrivialCone>())
          .def(CopyableVisitor<TrivialCone>());
      }
    };

  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_set_trivial_cones_hpp__
