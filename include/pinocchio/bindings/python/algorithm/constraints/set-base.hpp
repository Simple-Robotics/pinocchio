//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_set_base_hpp__
#define __pinocchio_python_algorithm_constraints_set_base_hpp__

#include <eigenpy/eigenpy.hpp>

#include "pinocchio/algorithm/constraints/box-set.hpp"

#include "pinocchio/bindings/python/utils/cast.hpp"
#include "pinocchio/bindings/python/utils/copyable.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    // Break ambiguity between SetBase::project(x) and SetBase::project(x, y)
    template<typename Set, typename VectorLike>
    static typename PINOCCHIO_EIGEN_PLAIN_TYPE(VectorLike)
      set_project_proxy(const Set & set, const Eigen::MatrixBase<VectorLike> & x)
    {
      return set.template project<VectorLike>(x);
    }

    template<typename Set, typename VectorLike>
    struct SetPythonVisitor : public boost::python::def_visitor<SetPythonVisitor<Set, VectorLike>>
    {
      typedef typename Set::Scalar Scalar;
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl.def(
            "project", set_project_proxy<Set, VectorLike>, bp::args("self", "f"),
            "Normal projection of a vector f onto the cone.")
          .def("dim", &Set::dim, "Returns the dimension of the cone.")
          .def("size", &Set::size, "Returns the size of the cone.")
#ifndef PINOCCHIO_PYTHON_SKIP_COMPARISON_OPERATIONS
          .def(
            "isInside",
            +[](const Set & self, const Eigen::MatrixBase<VectorLike> & f) -> bool {
              return self.template isInside<VectorLike>(f, Scalar(0));
            },
            bp::args("self", "f"), "Resize the constraint given active limits.")
          .def(bp::self == bp::self)
          .def(bp::self != bp::self)
#endif
          ;
      }
    };

    template<typename ConeSet>
    struct ConeSetPythonVisitor : public boost::python::def_visitor<ConeSetPythonVisitor<ConeSet>>
    {
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl.def("dual", &ConeSet::dual, bp::arg("self"), "Returns the dual cone associated to this");
      }
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_set_base_hpp__
