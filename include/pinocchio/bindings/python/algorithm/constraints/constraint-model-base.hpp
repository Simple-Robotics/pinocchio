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
    namespace internal
    {

      template<bool HasBaumgarteCorrector>
      struct BaumgarteCorrectorGetterSetter
      {
        template<typename VectorOut, typename VectorIn>
        static VectorOut getParameter(const VectorIn &)
        {
          assert(false && "Should not be reached");
          return ::pinocchio::visitors::internal::NoRun<VectorOut>::run();
        }

        template<typename VectorParameter, typename VectorNewParameter>
        static void setParameter(const VectorParameter &, const VectorNewParameter &)
        {
          assert(false && "Should not be reached");
        }
      };

      template<>
      struct BaumgarteCorrectorGetterSetter<true>
      {
        template<typename VectorOut, typename VectorIn>
        static VectorOut getParameter(const VectorIn & param)
        {
          return param;
        }

        template<typename VectorParameter, typename VectorNewParameter>
        static void setParameter(
          const Eigen::MatrixBase<VectorParameter> & param,
          const Eigen::MatrixBase<VectorNewParameter> & new_param)
        {
          param.const_cast_derived() = new_param;
        }
      };

    } // namespace internal

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
      typedef typename Self::ComplianceVectorTypeRef ComplianceVectorTypeRef;
      typedef typename Self::ComplianceVectorTypeConstRef ComplianceVectorTypeConstRef;

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
          .add_property(
            "baumgarte_corrector_parameters_kp",
            bp::make_function( //
              +[](const Self & self) -> context::VectorXs {
                static constexpr bool has_baumgarte_corrector =
                  traits<Self>::has_baumgarte_corrector;
                if (!has_baumgarte_corrector)
                {
                  std::stringstream ss;
                  ss << self.shortname() << " has no baumgarte corrector parameters.\n";
                  PINOCCHIO_THROW(std::invalid_argument, ss.str());
                }
                typedef internal::BaumgarteCorrectorGetterSetter<has_baumgarte_corrector> Getter;
                return Getter::template getParameter<context::VectorXs>(
                  self.baumgarte_corrector_parameters().Kp);
              }),
            bp::make_function( //
              +[](Self & self, const context::VectorXs & new_kp) {
                static constexpr bool has_baumgarte_corrector =
                  traits<Self>::has_baumgarte_corrector;
                if (!has_baumgarte_corrector)
                {
                  std::stringstream ss;
                  ss << self.shortname() << " has no baumgarte corrector parameters.\n";
                  PINOCCHIO_THROW(std::invalid_argument, ss.str());
                }
                typedef internal::BaumgarteCorrectorGetterSetter<has_baumgarte_corrector> Setter;
                return Setter::setParameter(self.baumgarte_corrector_parameters().Kp, new_kp);
              }),
            "Proportional term of baumgarte corrector of the constraint.")
          .add_property(
            "baumgarte_corrector_parameters_kd",
            bp::make_function( //
              +[](const Self & self) -> context::VectorXs {
                static constexpr bool has_baumgarte_corrector =
                  traits<Self>::has_baumgarte_corrector;
                if (!has_baumgarte_corrector)
                {
                  std::stringstream ss;
                  ss << self.shortname() << " has no baumgarte corrector parameters.\n";
                  PINOCCHIO_THROW(std::invalid_argument, ss.str());
                }
                typedef internal::BaumgarteCorrectorGetterSetter<has_baumgarte_corrector> Getter;
                return Getter::template getParameter<context::VectorXs>(
                  self.baumgarte_corrector_parameters().Kd);
              }),
            bp::make_function( //
              +[](Self & self, const context::VectorXs & new_kd) {
                static constexpr bool has_baumgarte_corrector =
                  traits<Self>::has_baumgarte_corrector;
                if (!has_baumgarte_corrector)
                {
                  std::stringstream ss;
                  ss << self.shortname() << " has no baumgarte corrector parameters.\n";
                  PINOCCHIO_THROW(std::invalid_argument, ss.str());
                }
                typedef internal::BaumgarteCorrectorGetterSetter<has_baumgarte_corrector> Setter;
                return Setter::setParameter(self.baumgarte_corrector_parameters().Kd, new_kd);
              }),
            "Damping term of baumgarte corrector of the constraint.")
          .def(
            "size", +[](const Self & self) -> int { return self.size(); }, "Constraint size.")
          .def(
            "calc", &calc, bp::args("self", "model", "data", "constraint_data"),
            "Evaluate the constraint values at the current state given by data and store the "
            "results.")
          .def(
            "jacobian", &jacobian, bp::args("self", "model", "data", "constraint_data"),
            "Compute the constraint jacobian.")
          .def(
            "jacobian_matrix_product", &jacobianMatrixProduct,
            bp::args("self", "model", "data", "constraint_data", "matrix"),
            "Forward chain rule: return product between the jacobian and a matrix.")
          .def(
            "jacobian_transpose_matrix_product", &jacobianTransposeMatrixProduct,
            bp::args("self", "model", "data", "constraint_data", "matrix"),
            "Backward chain rule: return product between the jacobian transpose and a matrix.")
          .def(
            "getRowSparsityPattern", &Self::getRowSparsityPattern, bp::args("self", "row_id"),
            bp::return_value_policy<bp::copy_const_reference>(),
            "Colwise sparsity associated with a given row.")
          .def(
            "getRowActiveIndexes", &Self::getRowActiveIndexes, bp::args("self", "row_id"),
            bp::return_value_policy<bp::copy_const_reference>(),
            "Vector of the active indexes associated with a given row.")
#ifndef PINOCCHIO_PYTHON_SKIP_COMPARISON_OPERATIONS
          .def(bp::self == bp::self)
          .def(bp::self != bp::self)
#endif
          ;
      }

      static void calc(
        const Self & self, const Model & model, const Data & data, ConstraintData & constraint_data)
      {
        self.calc(model, data, constraint_data);
      }

      static context::MatrixXs jacobian(
        const Self & self, const Model & model, const Data & data, ConstraintData & constraint_data)
      {
        const context::MatrixXs res(self.size(), model.nv);
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
        const context::MatrixXs res(self.size(), model.nv);
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
        const context::MatrixXs res(self.size(), model.nv);
        self.jacobianTransposeMatrixProduct(model, data, constraint_data, matrix, res);
        return res;
      }
    };
  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_constraints_model_base_hpp__
