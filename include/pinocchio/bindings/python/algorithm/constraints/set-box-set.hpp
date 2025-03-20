//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_set_box_set_hpp__
#define __pinocchio_python_algorithm_constraints_set_box_set_hpp__

#include <eigenpy/eigenpy.hpp>

#include "pinocchio/algorithm/constraints/box-set.hpp"

#include "pinocchio/bindings/python/algorithm/constraints/set-base.hpp"
#include "pinocchio/bindings/python/utils/cast.hpp"
#include "pinocchio/bindings/python/utils/copyable.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    template<typename BoxSet>
    struct BoxSetPythonVisitor : public boost::python::def_visitor<BoxSetPythonVisitor<BoxSet>>
    {
      typedef typename BoxSet::Scalar Scalar;
      typedef typename BoxSet::Vector Vector;
      typedef BoxSet Self;

      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl.def(bp::init<Eigen::DenseIndex>(
                 bp::args("self", "size"),
                 "Default constructor. By default, the bounds are set to ±inf."))
          .def(bp::init<const Self &>(bp::args("self", "other"), "Copy constructor."))
          .def(bp::init<context::VectorXs, context::VectorXs>(
            bp::args("self", "lb", "ub"), "Constructor from lower and upper bounds."))
          .def(
            "lb", (Vector & (Self::*)()) & Self::lb,
            "Returns a reference to the vector of lower bounds", bp::return_internal_reference<>())
          .def(
            "ub", (Vector & (Self::*)()) & Self::ub,
            "Returns a reference to the vector of upper bounds", bp::return_internal_reference<>())
          .def("resize", &BoxSet::resize, bp::args("self", "size"), "Resize the set.")
          .def(
            "conservativeResize", &BoxSet::conservativeResize, bp::args("self", "size"),
            "Resize the set following Eigen convention.")
          .def("isValid", &BoxSet::isValid, "Check if the constraint set is well defined.");
      }

      static void expose()
      {
        bp::class_<BoxSet>(
          "BoxSet", "Box set defined by a lower and an upper bounds [lb;ub].\n", bp::no_init)
          .def(SetPythonVisitor<BoxSet, context::VectorXs>())
          .def(BoxSetPythonVisitor())
          // .def(CastVisitor<BoxSet>())
          // .def(ExposeConstructorByCastVisitor<BoxSet,::pinocchio::BoxSet>())
          .def(CopyableVisitor<BoxSet>());
      }
    };

  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_set_box_set_hpp__
