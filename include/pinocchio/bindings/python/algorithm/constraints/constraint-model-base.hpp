//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_model_base_hpp__
#define __pinocchio_python_algorithm_constraints_model_base_hpp__

#include <boost/python.hpp>
#include <eigenpy/exception.hpp>
#include <eigenpy/eigen-to-python.hpp>

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-base.hpp"
#include "pinocchio/bindings/python/fwd.hpp"
#include "pinocchio/bindings/python/utils/macros.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    template<class ConstraintModelDerived>
    struct ConstraintModelBasePythonVisitor
    : public bp::def_visitor<ConstraintModelBasePythonVisitor<ConstraintModelDerived>>
    {
      typedef ConstraintModelDerived Self;
      typedef typename Self::Scalar Scalar;
      typedef typename Self::ConstraintSet ConstraintSet;
      typedef typename Self::ConstraintData ConstraintData;
      typedef context::Model Model;
      typedef context::Data Data;

    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl.PINOCCHIO_ADD_PROPERTY(Self, name, "Name of the constraint.")
          .def("classname", &Self::classname)
          .staticmethod("classname")
          .def("shortname", &Self::shortname, "Shortame for the constraint type")
          .def(
            "createData", &Self::createData, "Create a Data object for the given constraint model.")
          .add_property(
            "set",
            bp::make_function(
              (ConstraintSet & (Self::*)()) & Self::set, bp::return_internal_reference<>()),
            +[](Self & self, const ConstraintSet & new_set) { self.set() = new_set; },
            "Constraint set")
          .def(
            "size", +[](const Self & self) -> int { return self.size(); }, "Constraint size")
          .def("calc", &ConstraintModelBasePythonVisitor::_calc, bp::args("self", "model", "data", "constraint_data"))
          .def("jacobian", &ConstraintModelBasePythonVisitor::_jacobian, bp::args("self", "model", "data", "constraint_data"))
          .def("jacobian_matrix_product", &ConstraintModelBasePythonVisitor::_jacobianMatrixProduct,
            bp::args("self", "model", "data", "constraint_data", "matrix"))
          .def("jacobian_transpose_matrix_product", &ConstraintModelBasePythonVisitor::_jacobianTransposeMatrixProduct,
            bp::args("self", "model", "data", "constraint_data", "matrix"))
          .def(
            "compliance", +[](const Self & self) -> context::VectorXs { return self.compliance(); }, "Compliance")
          .def(
            "getRowSparsityPattern",
            &Self::getRowSparsityPattern,
            bp::args("self", "row_id"),
            bp::return_value_policy<bp::copy_const_reference>(),
            "Colwise sparsity associated with a given row.")
          .def(
            "getRowActiveIndexes",
            &Self::getRowActiveIndexes,
            bp::args("self", "row_id"),
            bp::return_value_policy<bp::copy_const_reference>(),
            "Vector of the active indexes associated with a given row.")
#ifndef PINOCCHIO_PYTHON_SKIP_COMPARISON_OPERATIONS
          .def(bp::self == bp::self)
          .def(bp::self != bp::self)
#endif
          ;
      }

      static void _calc(
        const Self & self, const Model & model, const Data & data, ConstraintData & constraint_data)
      {
        self.calc(model, data, constraint_data);
      }

      static context::MatrixXs _jacobian(
        const Self & self, const Model & model, const Data & data, ConstraintData & constraint_data)
      {
        const context::MatrixXs res(self.size(), model.nv);
        self.jacobian(model, data, constraint_data, res);
        return res;
      }

      static context::MatrixXs _jacobianMatrixProduct(
        const Self & self,
        const Model & model,
        const Data & data,
        const ConstraintData & constraint_data,
        const context::MatrixXs & matrix)
      {
        const context::MatrixXs res(self.size(), model.nv);
        self.jacobianMatrixProduct(model, data, constraint_data, matrix, res);
        return res;
      }

      static context::MatrixXs _jacobianTransposeMatrixProduct(
        const Self & self,
        const Model & model,
        const Data & data,
        const ConstraintData & constraint_data,
        const context::MatrixXs & matrix)
      {
        const context::MatrixXs res(self.size(), model.nv);
        self.jacobianTransposeMatrixProduct(model, data, constraint_data, matrix, res);
        return res;
      }
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_model_base_hpp__
