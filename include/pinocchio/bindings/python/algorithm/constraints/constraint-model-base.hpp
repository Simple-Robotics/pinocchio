//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_python_algorithm_constraints_model_base_hpp__
#define __pinocchio_python_algorithm_constraints_model_base_hpp__

#include <iostream>
#include <boost/python.hpp>
#include <eigenpy/exception.hpp>
#include <eigenpy/eigen-to-python.hpp>

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/algorithm/constraints/fwd.hpp"
#include "pinocchio/algorithm/constraints/constraint-model-base.hpp"

#include "pinocchio/bindings/python/fwd.hpp"
#include "pinocchio/bindings/python/utils/macros.hpp"
#include "pinocchio/bindings/python/utils/eigen.hpp"
#include "pinocchio/bindings/python/algorithm/constraints/baumgarte-corrector-vector-parameters.hpp"
#include "pinocchio/bindings/python/algorithm/constraints/baumgarte-corrector-parameters.hpp"

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
      typedef typename traits<Self>::ComplianceVectorTypeRef ComplianceVectorTypeRef;
      typedef typename traits<Self>::ComplianceVectorTypeConstRef ComplianceVectorTypeConstRef;

    public:
      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl.PINOCCHIO_ADD_PROPERTY(Self, name, "Name of the constraint.")
          .def("classname", &Self::classname)
          .staticmethod("classname")
          .def("shortname", &Self::shortname, "Short name of the class.")
          .def(
            "createData", &Self::createData, "Create a Data object for the given constraint model.")
          .add_property(
            "set",
            bp::make_function( //
              +[](const Self & self) -> const ConstraintSet & { return self.set(); },
              bp::return_internal_reference<>()),
            bp::make_function( //
              +[](Self & self, const ConstraintSet & new_set) { self.set() = new_set; }),
            "Constraint set.")
          .add_property(
            "compliance",
            bp::make_function( //
              +[](const Self & self) -> context::VectorXs { return self.compliance(); }),
            bp::make_function( //
              +[](Self & self, const context::VectorXs & new_vector) {
                self.compliance() = new_vector;
              }),
            "Compliance of the constraint.")
          .def(
            "size", +[](const Self & self) -> int { return self.size(); }, "Constraint size.")
          .def(
            "activeSize", +[](const Self & self) -> int { return self.activeSize(); },
            "Constraint active size.")
          .def(
            "calc", &calc, bp::args("self", "model", "data", "constraint_data"),
            "Evaluate the constraint values at the current state given by data and store the "
            "results.")
          .def(
            "jacobian", &jacobian, bp::args("self", "model", "data", "constraint_data"),
            "Compute the constraint jacobian.")
          .def(
            "jacobianMatrixProduct", &jacobianMatrixProduct,
            bp::args("self", "model", "data", "constraint_data", "matrix"),
            "Forward chain rule: return product between the jacobian and a matrix.")
          .def(
            "jacobianTransposeMatrixProduct", &jacobianTransposeMatrixProduct,
            bp::args("self", "model", "data", "constraint_data", "matrix"),
            "Backward chain rule: return product between the jacobian transpose and a matrix.")
          .def(
            "getRowActivableSparsityPattern", &Self::getRowActivableSparsityPattern,
            bp::args("self", "row_id"), bp::return_value_policy<bp::copy_const_reference>(),
            "Colwise sparsity associated with a given row.")
          .def(
            "getRowActiveSparsityPattern", &Self::getRowActiveSparsityPattern,
            bp::args("self", "row_id"), bp::return_value_policy<bp::copy_const_reference>(),
            "Active colwise sparsity associated with a given row.")
          .def(
            "getRowActivableIndexes", &Self::getRowActivableIndexes, bp::args("self", "row_id"),
            bp::return_value_policy<bp::copy_const_reference>(),
            "Vector of the activable indexes associated with a given row.")
          .def(
            "getRowActiveIndexes", &Self::getRowActiveIndexes, bp::args("self", "row_id"),
            bp::return_value_policy<bp::copy_const_reference>(),
            "Vector of the active indexes associated with a given row.")
          .def(
            "getActiveCompliance", bp::make_function(+[](const Self & self) -> context::VectorXs {
              return self.getActiveCompliance();
            }),
            "Vector of the active compliance internally stored in the constraint.")
#ifndef PINOCCHIO_PYTHON_SKIP_COMPARISON_OPERATIONS
          .def(bp::self == bp::self)
          .def(bp::self != bp::self)
#endif
          ;
        // CHOICE: right now we use the scalar Baumgarte
        // if (::pinocchio::traits<ConstraintModelDerived>::has_baumgarte_corrector_vector)
        // {
        //   typedef typename traits<ConstraintModelDerived>::BaumgarteCorrectorVectorParameters
        //     BaumgarteCorrectorVectorParameters;
        //   typedef typename traits<ConstraintModelDerived>::BaumgarteCorrectorVectorParametersRef
        //     BaumgarteCorrectorVectorParametersRef;

        //   typedef typename std::conditional<
        //     std::is_reference<BaumgarteCorrectorVectorParametersRef>::value,
        //     bp::return_internal_reference<>, bp::with_custodian_and_ward_postcall<0, 1>>::type
        //     ReturnPolicy;

        //   cl.add_property(
        //     "baumgarte_corrector_vector_parameters",
        //     bp::make_function( //
        //       +[](Self & self) -> BaumgarteCorrectorVectorParametersRef {
        //         return self.baumgarte_corrector_vector_parameters();
        //       },
        //       ReturnPolicy()),
        //     bp::make_function( //
        //       +[](Self & self, const BaumgarteCorrectorVectorParameters & copy) {
        //         self.baumgarte_corrector_vector_parameters() = copy;
        //       }),
        //     "Baumgarte vector parameters associated with the constraint.");

        //   typedef typename BaumgarteCorrectorVectorParameters::VectorType BaumgarteVectorType;
        //   const std::string BaumgarteVectorType_name = getEigenTypeName<BaumgarteVectorType>();

        //   const std::string BaumgarteCorrectorVectorParameter_classname =
        //     "BaumgarteCorrectorVectorParameters_" + BaumgarteVectorType_name;

        //   BaumgarteCorrectorVectorParametersPythonVisitor<BaumgarteCorrectorVectorParameters>::
        //     expose(BaumgarteCorrectorVectorParameter_classname);
        // }
        if (::pinocchio::traits<ConstraintModelDerived>::has_baumgarte_corrector)
        {
          typedef BaumgarteCorrectorParametersTpl<Scalar> BaumgarteCorrectorParameters;
          BaumgarteCorrectorParametersPythonVisitor<BaumgarteCorrectorParameters>::expose();

          cl.add_property(
            "baumgarte_corrector_parameters",
            bp::make_function( //
              +[](Self & self) -> BaumgarteCorrectorParameters & {
                return self.baumgarte_corrector_parameters();
              },
              bp::return_internal_reference<>()),
            bp::make_function( //
              +[](Self & self, const BaumgarteCorrectorParameters & copy) {
                self.baumgarte_corrector_parameters() = copy;
              },
              bp::return_internal_reference<>()),
            "Baumgarte parameters associated with the constraint.");
        }
      }

      static void calc(
        const Self & self, const Model & model, const Data & data, ConstraintData & constraint_data)
      {
        self.calc(model, data, constraint_data);
      }

      static context::MatrixXs jacobian(
        const Self & self, const Model & model, const Data & data, ConstraintData & constraint_data)
      {
        context::MatrixXs res = context::MatrixXs::Zero(self.size(), model.nv);
        self.jacobian(model, data, constraint_data, res);
        return res;
      }

      static context::MatrixXs jacobianMatrixProduct(
        const Self & self,
        const Model & model,
        const Data & data,
        const ConstraintData & constraint_data,
        const context::MatrixXs & matrix)
      {
        context::MatrixXs res = context::MatrixXs::Zero(self.size(), model.nv);
        self.jacobianMatrixProduct(model, data, constraint_data, matrix, res);
        return res;
      }

      static context::MatrixXs jacobianTransposeMatrixProduct(
        const Self & self,
        const Model & model,
        const Data & data,
        const ConstraintData & constraint_data,
        const context::MatrixXs & matrix)
      {
        context::MatrixXs res = context::MatrixXs::Zero(self.size(), model.nv);
        self.jacobianTransposeMatrixProduct(model, data, constraint_data, matrix, res);
        return res;
      }
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_model_base_hpp__
