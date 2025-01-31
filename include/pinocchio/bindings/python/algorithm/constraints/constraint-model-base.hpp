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
      typedef ConstraintModelDerived Self;
      typedef typename ConstraintModelDerived::Scalar Scalar;

      typedef ModelTpl<Scalar, ConstraintModelDerived::Options, JointCollectionDefaultTpl> Model;
      typedef DataTpl<Scalar, ConstraintModelDerived::Options, JointCollectionDefaultTpl> Data;

      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl
          .PINOCCHIO_ADD_PROPERTY(Self, name, "Name of the contact.")
          .def("classname", &Self::classname)
          .staticmethod("classname")
          .def("createData", &Self::createData, bp::arg("self") ,
            "Create a Data object for the given constraint model.")
          .def("shortname", &Self::shortname, bp::arg("self") ,
            "Shortame for the constraint type")
          .def("calc", &calc, bp::arg("self", "model", "data", "constraint_data"))
          .def("jacobian", &jacobian, bp::arg("self", "model", "data", "jacobian_matrix"))
          // .def("jacobianMatrixProduct", ...)
          // .def("jacobianTransposeMatrixProduct", ...)
#ifndef PINOCCHIO_PYTHON_SKIP_COMPARISON_OPERATIONS
          .def(bp::self == bp::self)
          .def(bp::self != bp::self)
#endif
          ;
      }

      static void calc(const Self & self, Model & model, Data & data)
      {
        return self.calc(model, data, constraint_data);
      }

      static void jacobian(
        const Self & self, Model & model, Data & data, ConstraintData & constraint_data,
        context::MatrixXs & jacobian_matrix)
      {
        return self.jacobian(model, data, constraint_data, jacobian_matrix);
      }
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_model_base_hpp__
