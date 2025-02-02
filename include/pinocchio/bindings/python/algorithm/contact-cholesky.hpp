//
// Copyright (c) 2020-2025 INRIA
//

#ifndef __pinocchio_python_algorithm_contact_cholesky_hpp__
#define __pinocchio_python_algorithm_contact_cholesky_hpp__

#include <eigenpy/memory.hpp>
#include "pinocchio/algorithm/contact-cholesky.hpp"

#include "pinocchio/algorithm/delassus-operator-dense.hpp"
#include "pinocchio/algorithm/delassus-operator-sparse.hpp"

#include "pinocchio/bindings/python/utils/macros.hpp"
#include "pinocchio/bindings/python/utils/std-vector.hpp"
#include "pinocchio/bindings/python/utils/comparable.hpp"
#include "pinocchio/bindings/python/utils/copyable.hpp"

#include "pinocchio/bindings/python/algorithm/delassus-operator.hpp"

namespace pinocchio
{
  namespace python
  {
    namespace bp = boost::python;

    template<typename ContactCholeskyDecomposition>
    struct ContactCholeskyDecompositionPythonVisitor
    : public boost::python::def_visitor<
        ContactCholeskyDecompositionPythonVisitor<ContactCholeskyDecomposition>>
    {
      typedef ContactCholeskyDecomposition Self;
      typedef typename Self::Scalar Scalar;
      typedef typename Self::RigidConstraintModel RigidConstraintModel;
      typedef typename Self::RigidConstraintData RigidConstraintData;
      typedef typename Self::Matrix Matrix;
      typedef typename Self::RowMatrix RowMatrix;
      typedef typename Self::Vector Vector;

      typedef context::ConstraintModel ConstraintModel;
      typedef context::ConstraintData ConstraintData;

      typedef Eigen::Ref<Matrix> RefMatrix;
      typedef Eigen::Ref<RowMatrix> RefRowMatrix;
      typedef Eigen::Ref<Vector> RefVector;

      typedef typename PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintModel)
        RigidConstraintModelVector;
      typedef typename PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(RigidConstraintData)
        RigidConstraintDataVector;

      typedef
        typename PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(ConstraintModel) ConstraintModelVector;
      typedef
        typename PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(ConstraintData) ConstraintDataVector;

      typedef pinocchio::python::context::Model Model;
      typedef pinocchio::python::context::Data Data;

      template<class PyClass>
      void visit(PyClass & cl) const
      {
        cl.def(bp::init<>(bp::arg("self"), "Default constructor."))
          .def(bp::init<const Model &>(bp::args("self", "model"), "Constructor from a model."))
          .def(bp::init<const Model &, const RigidConstraintModelVector &>(
            (bp::arg("self"), bp::arg("model"), bp::arg("contact_models")),
            "Constructor from a model and a collection of RigidConstraintModels."))
          .def(bp::init<const Model &, const ConstraintModelVector &>(
            (bp::arg("self"), bp::arg("model"), bp::arg("constraint_models")),
            "Constructor from a model and a collection of ConstraintModels."))

          .add_property(
            "U",
            bp::make_function(
              +[](const Self & self) -> RefRowMatrix { return RefRowMatrix(self.U); },
              bp::with_custodian_and_ward_postcall<0, 1>()),
            "")
          .add_property(
            "D",
            bp::make_function(
              +[](const Self & self) -> RefVector { return RefVector(self.D); },
              bp::with_custodian_and_ward_postcall<0, 1>()),
            "")
          .add_property(
            "Dinv",
            bp::make_function(
              +[](const Self & self) -> RefVector { return RefVector(self.Dinv); },
              bp::with_custodian_and_ward_postcall<0, 1>()),
            "")

          .def("size", &Self::size, bp::arg("self"), "Size of the decomposition.")
          .def(
            "constraintDim", &Self::constraintDim, bp::arg("self"),
            "Returns the total dimension of the constraints contained in the Cholesky "
            "factorization.")
          .def(
            "matrix", (Matrix(Self::*)(void) const)&Self::matrix, bp::arg("self"),
            "Returns the matrix resulting from the decomposition.")

          .def(
            "compute",
            (void (*)(
              Self & self, const Model &, Data &, const RigidConstraintModelVector &,
              RigidConstraintDataVector &, const Scalar))&compute,
            (bp::arg("self"), bp::arg("model"), bp::arg("data"), bp::arg("contact_models"),
             bp::arg("contact_datas"), bp::arg("mu") = 0),
            "Computes the Cholesky decompostion of the augmented matrix containing the KKT matrix\n"
            "related to the system mass matrix and the Jacobians of the contact patches contained "
            "in\n"
            "the vector of RigidConstraintModel named contact_models. The decomposition is "
            "regularized with a factor mu.")

          .def(
            "compute",
            (void (*)(
              Self & self, const Model &, Data &, const RigidConstraintModelVector &,
              RigidConstraintDataVector &, const Vector &))&compute,
            (bp::arg("self"), bp::arg("model"), bp::arg("data"), bp::arg("contact_models"),
             bp::arg("contact_datas"), bp::arg("mus")),
            "Computes the Cholesky decompostion of the augmented matrix containing the KKT matrix\n"
            "related to the system mass matrix and the Jacobians of the contact patches contained "
            "in\n"
            "the vector of RigidConstraintModel named contact_models. The decomposition is "
            "regularized with a factor mu.")

          .def(
            "compute",
            (void (*)(
              Self & self, const Model &, Data &, const ConstraintModelVector &,
              ConstraintDataVector &, const Scalar))&compute,
            (bp::arg("self"), bp::arg("model"), bp::arg("data"), bp::arg("constraint_models"),
             bp::arg("constraint_datas"), bp::arg("mu") = 0),
            "Computes the Cholesky decompostion of the augmented matrix containing the KKT matrix\n"
            "related to the system mass matrix and the Jacobians of the contact patches contained "
            "in\n"
            "the vector of ConstraintModel named constraint_models. The decomposition is "
            "regularized with a factor mu.")

          .def(
            "compute",
            (void (*)(
              Self & self, const Model &, Data &, const ConstraintModelVector &,
              ConstraintDataVector &, const Vector &))&compute,
            (bp::arg("self"), bp::arg("model"), bp::arg("data"), bp::arg("constraint_models"),
             bp::arg("constraint_datas"), bp::arg("mus")),
            "Computes the Cholesky decompostion of the augmented matrix containing the KKT matrix\n"
            "related to the system mass matrix and the Jacobians of the contact patches contained "
            "in\n"
            "the vector of ConstraintModel named constraint_models. The decomposition is "
            "regularized with a factor mu.")

          .def(
            "updateDamping", (void(Self::*)(const Scalar &)) & Self::updateDamping,
            bp::args("self", "mu"),
            "Update the damping term on the upper left block part of the KKT matrix. The damping "
            "term should be positive.")

          .def(
            "updateDamping", &Self::template updateDamping<Vector>, bp::args("self", "mus"),
            "Update the damping terms on the upper left block part of the KKT matrix. The damping "
            "terms should be all positives.")

          .def(
            "getDamping", +[](const Self & self) -> Vector { return Vector(self.getDamping()); },
            bp::arg("self"), "Returns the current damping vector.")

          .def(
            "getInverseOperationalSpaceInertiaMatrix",
            (Matrix(Self::*)(bool) const)&Self::getInverseOperationalSpaceInertiaMatrix,
            (bp::arg("self"), bp::arg("enforce_symmetry") = false),
            "Returns the Inverse of the Operational Space Inertia Matrix resulting from the "
            "decomposition.",
            bp::return_value_policy<bp::return_by_value>())

          .def(
            "getOperationalSpaceInertiaMatrix",
            (Matrix(Self::*)(void) const)&Self::getOperationalSpaceInertiaMatrix, bp::arg("self"),
            "Returns the Operational Space Inertia Matrix resulting from the decomposition.",
            bp::return_value_policy<bp::return_by_value>())

          .def(
            "getInverseMassMatrix", (Matrix(Self::*)(void) const)&Self::getInverseMassMatrix,
            bp::arg("self"),
            "Returns the inverse of the Joint Space Inertia Matrix or \"mass matrix\".",
            bp::return_value_policy<bp::return_by_value>())

          .def(
            "solve", &solve<Matrix>, bp::args("self", "matrix"),
            "Computes the solution of \f$ A x = b \f$ where self corresponds to the Cholesky "
            "decomposition of A.",
            bp::return_value_policy<bp::return_by_value>())

          .def(
            "inverse", (Matrix(Self::*)(void) const)&Self::inverse, bp::arg("self"),
            "Returns the inverse matrix resulting from the decomposition.")

          .def(
            "getMassMatrixChoeslkyDecomposition",
            &Self::template getMassMatrixChoeslkyDecomposition<
              Scalar, 0, JointCollectionDefaultTpl>,
            bp::arg("self"),
            "Retrieves the Cholesky decomposition of the Mass Matrix contained in the current "
            "decomposition.")

          .def(
            "getDelassusCholeskyExpression", &Self::getDelassusCholeskyExpression, bp::arg("self"),
            "Returns the Cholesky decomposition expression associated to the underlying Delassus "
            "matrix.",
            bp::with_custodian_and_ward_postcall<0, 1>())

          .def(ComparableVisitor<Self, pinocchio::is_floating_point<Scalar>::value>());
      }

      static void expose()
      {
        bp::class_<ContactCholeskyDecomposition>(
          "ContactCholeskyDecomposition",
          "Contact information container for contact dynamic algorithms.", bp::no_init)
          .def(ContactCholeskyDecompositionPythonVisitor<ContactCholeskyDecomposition>())
          .def(CopyableVisitor<ContactCholeskyDecomposition>());

        {
          typedef typename ContactCholeskyDecomposition::DelassusCholeskyExpression
            DelassusCholeskyExpression;

          bp::class_<DelassusCholeskyExpression>(
            "DelassusCholeskyExpression",
            "Delassus Cholesky expression associated to a given ContactCholeskyDecomposition "
            "object.",
            bp::no_init)
            .def(bp::init<const ContactCholeskyDecomposition &>(
              bp::args("self", "cholesky_decomposition"),
              "Build from a given ContactCholeskyDecomposition object.")
                   [bp::with_custodian_and_ward<1, 2>()])
            .def(
              "cholesky",
              +[](const DelassusCholeskyExpression & self) -> ContactCholeskyDecomposition & {
                return const_cast<ContactCholeskyDecomposition &>(self.cholesky());
              },
              bp::arg("self"),
              "Returns the Constraint Cholesky decomposition associated to this "
              "DelassusCholeskyExpression.",
              bp::return_internal_reference<>())

            .def(DelassusOperatorBasePythonVisitor<DelassusCholeskyExpression>());
        }

        {
          typedef DelassusOperatorDenseTpl<context::Scalar, context::Options> DelassusOperatorDense;
          bp::class_<DelassusOperatorDense>(
            "DelassusOperatorDense", "Delassus Cholesky dense operator from a dense matrix.",
            bp::no_init)
            .def(bp::init<const Eigen::Ref<const context::MatrixXs>>(
              bp::args("self", "matrix"), "Build from a given dense matrix"))

            .def(DelassusOperatorBasePythonVisitor<DelassusOperatorDense>());
        }

        {
          typedef DelassusOperatorSparseTpl<context::Scalar, context::Options>
            DelassusOperatorSparse;
          bp::class_<DelassusOperatorSparse, boost::noncopyable>(
            "DelassusOperatorSparse", "Delassus Cholesky sparse operator from a sparse matrix.",
            bp::no_init)
            .def(bp::init<const context::SparseMatrix &>(
              bp::args("self", "matrix"), "Build from a given sparse matrix"))

            .def(DelassusOperatorBasePythonVisitor<DelassusOperatorSparse>());
        }
#ifdef PINOCCHIO_WITH_ACCELERATE_SUPPORT
        {
          typedef Eigen::AccelerateLLT<context::SparseMatrix> AccelerateCholeskyDecomposition;
          typedef DelassusOperatorSparseTpl<
            context::Scalar, context::Options, AccelerateCholeskyDecomposition>
            DelassusOperatorSparseAccelerate;
          bp::class_<DelassusOperatorSparseAccelerate, boost::noncopyable>(
            "DelassusOperatorSparseAccelerate",
            "Delassus Cholesky sparse operator leveraging the Accelerate framework on APPLE "
            "systems.",
            bp::no_init)
            .def(bp::init<const context::SparseMatrix &>(
              bp::args("self", "matrix"), "Build from a given sparse matrix"))

            .def(DelassusOperatorBasePythonVisitor<DelassusOperatorSparseAccelerate>());
        }
#endif
      }

      template<typename MatrixType>
      static Matrix solve(const Self & self, const MatrixType & mat)
      {
        return self.solve(mat);
      }

      template<
        typename ConstraintModel,
        typename ConstraintModelAllocator,
        typename ConstraintData,
        typename ConstraintDataAllocator>
      static void compute(
        Self & self,
        const Model & model,
        Data & data,
        const std::vector<ConstraintModel, ConstraintModelAllocator> & contact_models,
        std::vector<ConstraintData, ConstraintDataAllocator> & contact_datas,
        const Scalar mu = static_cast<Scalar>(0))
      {
        self.compute(model, data, contact_models, contact_datas, mu);
      }

      template<
        typename ConstraintModel,
        typename ConstraintModelAllocator,
        typename ConstraintData,
        typename ConstraintDataAllocator>
      static void compute(
        Self & self,
        const Model & model,
        Data & data,
        const std::vector<ConstraintModel, ConstraintModelAllocator> & contact_models,
        std::vector<ConstraintData, ConstraintDataAllocator> & contact_datas,
        const Vector & mus)
      {
        self.compute(model, data, contact_models, contact_datas, mus);
      }
    };

  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_algorithm_contact_cholesky_hpp__
