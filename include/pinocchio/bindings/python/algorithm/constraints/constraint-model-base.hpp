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
      typedef typename ConstraintModelDerived::Scalar Scalar;
      typedef ModelTpl<Scalar, ConstraintModelDerived::Options, JointCollectionDefaultTpl> Model;
      typedef DataTpl<Scalar, ConstraintModelDerived::Options, JointCollectionDefaultTpl> Data;

    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl.PINOCCHIO_ADD_PROPERTY(Self, name, "Name of the contact.")
          .def("classname", &Self::classname)
          .staticmethod("classname")
          .def(
            "createData", &Self::createData, "Create a Data object for the given constraint model.")
          .def("shortname", &Self::shortname, "Shortame for the constraint type")
          .def("set", &Self::shortname, "Constraint set")
          .def("size", &Self::size, "Constraint size")
        // .def("compliance", &Self::compliance,
        //   "Return the compliance stored in the model.")
        // .def("calc", &calc, bp::args("self", "model", "data", "constraint_data"))
        // .def("jacobian", &jacobian, bp::args("self", "model", "data", "jacobian_matrix"))
        // .def("jacobian_matrix_product", &jacobianMatrixProduct,
        //   bp::args("self", "model", "data", "matrix"))
        // .def("jacobian_transpose_matrix_product", &jacobianTransposeMatrixProduct,
        //   bp::args("self", "model", "data", "matrix"))
#ifndef PINOCCHIO_PYTHON_SKIP_COMPARISON_OPERATIONS
          .def(bp::self == bp::self)
          .def(bp::self != bp::self)
#endif
          ;
      }

      static void calc(
        const Self & self, const Model & model, const Data & data, ConstraintData & constraint_data)
      {
        return self.calc(model, data, constraint_data);
      }

      static context::MatrixXs jacobian(
        const Self & self, const Model & model, const Data & data, ConstraintData & constraint_data)
      {
        const context::MatrixXs res(self.size(), model.nv);
        self.jacobian(model, data, constraint_data, res);
        return res;
      }

      static void jacobianMatrixProduct(
        const Self & self,
        Model & model,
        Data & data,
        ConstraintData & constraint_data,
        context::MatrixXs & matrix)
      {
        return self.jacobianMatrixProduct(model, data, constraint_data, matrix);
      }

      static void jacobianTransposeMatrixProduct(
        const Self & self,
        Model & model,
        Data & data,
        ConstraintData & constraint_data,
        context::MatrixXs & matrix)
      {
        return self.jacobianTransposeMatrixProduct(model, data, constraint_data, matrix);
      }
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_model_base_hpp__
