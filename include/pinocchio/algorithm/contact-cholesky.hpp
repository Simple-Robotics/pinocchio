//
// Copyright (c) 2019-2025 INRIA
//

#ifndef __pinocchio_algorithm_contact_cholesky_hpp__
#define __pinocchio_algorithm_contact_cholesky_hpp__

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/math/matrix-block.hpp"
#include "pinocchio/math/triangular-matrix.hpp"
#include "pinocchio/container/storage.hpp"

#include "pinocchio/algorithm/constraints/constraints.hpp"
#include <functional>

namespace pinocchio
{

  // Forward declaration of algo
  namespace details
  {
    template<typename MatrixLike, int ColsAtCompileTime = MatrixLike::ColsAtCompileTime>
    struct UvAlgo;

    template<typename MatrixLike, int ColsAtCompileTime = MatrixLike::ColsAtCompileTime>
    struct UtvAlgo;

    template<typename MatrixLike, int ColsAtCompileTime = MatrixLike::ColsAtCompileTime>
    struct UivAlgo;

    template<typename MatrixLike, int ColsAtCompileTime = MatrixLike::ColsAtCompileTime>
    struct UtivAlgo;

    template<typename Scalar, int Options, typename VectorLike>
    VectorLike & inverseAlgo(
      const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
      const Eigen::DenseIndex col,
      const Eigen::MatrixBase<VectorLike> & vec);
  } // namespace details

  template<typename _ContactCholeskyDecomposition>
  struct DelassusCholeskyExpressionTpl;

  ///
  ///  \brief Contact Cholesky decomposition structure. This structure allows
  ///        to compute in a efficient and parsimonious way the Cholesky decomposition
  ///        of the KKT matrix related to the contact dynamics.
  ///        Such a decomposition is usefull when computing both the forward dynamics in contact
  ///        or the related analytical derivatives.
  ///
  ///
  /// \tparam _Scalar Scalar type.
  ///  \tparam _Options Alignment Options of the Eigen objects contained in the data structure.
  ///
  template<typename _Scalar, int _Options>
  struct PINOCCHIO_UNSUPPORTED_MESSAGE("The API will change towards more flexibility")
    ContactCholeskyDecompositionTpl
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef pinocchio::Index Index;
    typedef _Scalar Scalar;
    enum
    {
      LINEAR = 0,
      ANGULAR = 3,
      Options = _Options
    };

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> Vector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> Matrix;
    typedef typename PINOCCHIO_EIGEN_PLAIN_ROW_MAJOR_TYPE(Matrix) RowMatrix;

    typedef EigenStorageTpl<Vector> EigenStorageVector;
    typedef EigenStorageTpl<Matrix> EigenStorageMatrix;
    typedef EigenStorageTpl<RowMatrix> EigenStorageRowMatrix;

    // TODO Remove when API is stabilized
    PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH
    PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS
    typedef RigidConstraintModelTpl<Scalar, Options> RigidConstraintModel;
    typedef RigidConstraintDataTpl<Scalar, Options> RigidConstraintData;
    PINOCCHIO_COMPILER_DIAGNOSTIC_POP
    typedef Eigen::Matrix<Eigen::DenseIndex, Eigen::Dynamic, 1, Options> EigenIndexVector;
    typedef
      typename PINOCCHIO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(EigenIndexVector) VectorOfEigenIndexVector;
    typedef Eigen::Matrix<bool, Eigen::Dynamic, 1, Options> BooleanVector;

    ///@{
    /// \brief Data information related to the Sparsity structure of the Cholesky decompostion
    struct Slice
    {
      Slice(const Eigen::DenseIndex & first_index, const Eigen::DenseIndex & size)
      : first_index(first_index)
      , size(size)
      {
      }

      Eigen::DenseIndex first_index;
      Eigen::DenseIndex size;
    };

    typedef DelassusCholeskyExpressionTpl<ContactCholeskyDecompositionTpl>
      DelassusCholeskyExpression;
    friend struct DelassusCholeskyExpressionTpl<ContactCholeskyDecompositionTpl>;

    typedef std::vector<Slice> SliceVector;
    typedef std::vector<SliceVector> VectorOfSliceVector;
    ///@}

    ///
    /// \brief Default constructor
    ///
    ContactCholeskyDecompositionTpl()
    : D(D_storage.map())
    , Dinv(Dinv_storage.map())
    , U(U_storage.map())
    , DUt(DUt_storage.map())
    , compliance(compliance_storage.map())
    , damping(damping_storage.map())
    {
    }

    ///
    /// \brief Constructor from a model.
    ///
    /// \param[in] model Model of the kinematic tree.
    ///
    template<typename S1, int O1, template<typename, int> class JointCollectionTpl>
    explicit ContactCholeskyDecompositionTpl(const ModelTpl<S1, O1, JointCollectionTpl> & model)
    : D(D_storage.map())
    , Dinv(Dinv_storage.map())
    , U(U_storage.map())
    , DUt(DUt_storage.map())
    , compliance(compliance_storage.map())
    , damping(damping_storage.map())
    {
      // TODO Remove when API is stabilized
      PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH
      PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS
      std::vector<ConstraintModel> empty_contact_models;
      PINOCCHIO_COMPILER_DIAGNOSTIC_POP
      resize(model, empty_contact_models);
    }

    ///
    /// \brief Constructor from a model and a collection of RigidConstraintModel objects.
    ///
    /// \param[in] model Model of the kinematic tree
    /// \param[in] constraint_models Vector of ConstraintModels
    /// information
    ///
    template<
      typename S1,
      int O1,
      template<typename, int> class JointCollectionTpl,
      class ConstraintModel,
      class ConstraintAllocator>
    ContactCholeskyDecompositionTpl(
      const ModelTpl<S1, O1, JointCollectionTpl> & model,
      const std::vector<ConstraintModel, ConstraintAllocator> & contact_models)
    : D(D_storage.map())
    , Dinv(Dinv_storage.map())
    , U(U_storage.map())
    , DUt(DUt_storage.map())
    , compliance(compliance_storage.map())
    , damping(damping_storage.map())
    {
      typedef std::reference_wrapper<const ConstraintModel> WrappedType;
      typedef std::vector<WrappedType> WrappedTypeVector;

      WrappedTypeVector wrapped_contact_models(contact_models.cbegin(), contact_models.cend());
      resize(model, wrapped_contact_models);
    }

    ///
    /// \brief Constructor from a model and a collection of RigidConstraintModel objects.
    ///
    /// \param[in] model Model of the kinematic tree
    /// \param[in] contact_models Vector of RigidConstraintModel references containing the contact
    /// information
    ///
    template<
      typename S1,
      int O1,
      template<typename, int> class JointCollectionTpl,
      template<typename T> class Holder,
      class Allocator>
    ContactCholeskyDecompositionTpl(
      const ModelTpl<S1, O1, JointCollectionTpl> & model,
      const std::vector<Holder<const RigidConstraintModelTpl<S1, O1>>, Allocator> & contact_models)
    : D(D_storage.map())
    , Dinv(Dinv_storage.map())
    , U(U_storage.map())
    , DUt(DUt_storage.map())
    , compliance(compliance_storage.map())
    , damping(damping_storage.map())
    {
      resize(model, contact_models);
    }

    ///
    /// \brief Copy constructor
    ///
    /// \param[in] other ContactCholeskyDecompositionTpl to copy
    ///
    ContactCholeskyDecompositionTpl(const ContactCholeskyDecompositionTpl & other)
    : D(D_storage.map())
    , Dinv(Dinv_storage.map())
    , U(U_storage.map())
    , DUt(DUt_storage.map())
    , compliance(compliance_storage.map())
    , damping(damping_storage.map())
    {
      *this = other;
    }

    ContactCholeskyDecompositionTpl & operator=(const ContactCholeskyDecompositionTpl & other)
    {
      parents_fromRow = other.parents_fromRow;
      nv_subtree_fromRow = other.nv_subtree_fromRow;
      nv = other.nv;

      rowise_sparsity_pattern = other.rowise_sparsity_pattern;

      D_storage = other.D_storage;
      Dinv_storage = other.Dinv_storage;
      U_storage = other.U_storage;
      DUt_storage = other.DUt_storage;
      compliance_storage = other.compliance_storage;
      damping_storage = other.damping_storage;

      return *this;
    }

    ///
    ///  \brief Internal memory allocation.
    ///
    /// \param[in] model Model of the kinematic tree
    /// \param[in] contact_models Vector of ConstraintModel
    ///
    template<
      typename S1,
      int O1,
      template<typename, int> class JointCollectionTpl,
      class ConstraintModel,
      class ConstraintAllocator>
    void resize(
      const ModelTpl<S1, O1, JointCollectionTpl> & model,
      const std::vector<ConstraintModel, ConstraintAllocator> & contact_models)
    {
      typedef std::reference_wrapper<const ConstraintModel> WrappedConstraintModelType;
      typedef std::vector<WrappedConstraintModelType> WrappedConstraintModelTypeVector;

      WrappedConstraintModelTypeVector wrapped_contact_models(
        contact_models.cbegin(), contact_models.cend());

      resize(model, wrapped_contact_models);
    }

    ///
    ///  \brief Internal memory allocation.
    ///
    /// \param[in] model Model of the kinematic tree
    /// \param[in] contact_models Vector of ConstraintModel
    /// \param[in] contact_datas Vector of ConstraintData
    ///
    template<
      typename S1,
      int O1,
      template<typename, int> class JointCollectionTpl,
      class ConstraintModel,
      class ConstraintAllocator,
      class ConstraintData,
      class ConstraintDataAllocator>
    PINOCCHIO_DEPRECATED void allocate(
      const ModelTpl<S1, O1, JointCollectionTpl> & model,
      const std::vector<ConstraintModel, ConstraintAllocator> & contact_models)
    {
      resize(model, contact_models);
    }

    ///
    ///  \brief Internal memory allocation.
    ///
    /// \param[in] model Model of the kinematic tree
    /// \param[in] contact_models Vector of ConstraintModel
    ///
    template<
      typename S1,
      int O1,
      template<typename, int> class JointCollectionTpl,
      template<typename T> class Holder,
      class ConstraintModel,
      class ConstraintAllocator>
    void resize(
      const ModelTpl<S1, O1, JointCollectionTpl> & model,
      const std::vector<Holder<const ConstraintModel>, ConstraintAllocator> & contact_models);

    ///
    /// \brief Returns the Inverse of the Operational Space Inertia Matrix resulting from the
    /// decomposition.
    ///
    Matrix getInverseOperationalSpaceInertiaMatrix(bool enforce_symmetry = false) const
    {
      Matrix res(constraintDim(), constraintDim());
      getInverseOperationalSpaceInertiaMatrix(res, enforce_symmetry);
      return res;
    }

    template<typename MatrixType>
    void getInverseOperationalSpaceInertiaMatrix(
      const Eigen::MatrixBase<MatrixType> & res, bool enforce_symmetry = false) const
    {
      const auto U1 = U.topLeftCorner(constraintDim(), constraintDim());

      const auto dim = constraintDim();
      typedef Eigen::Map<RowMatrix> MapRowMatrix;
      MapRowMatrix OSIMinv = MapRowMatrix(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, dim, dim));

      PINOCCHIO_EIGEN_MALLOC_NOT_ALLOWED();
      MatrixType & res_ = res.const_cast_derived();
      OSIMinv.noalias() = D.head(dim).asDiagonal() * U1.adjoint();
      res_.noalias() = -U1 * OSIMinv;
      if (enforce_symmetry)
        enforceSymmetry(res_);
      PINOCCHIO_EIGEN_MALLOC_ALLOWED();
    }

    /// \brief Returns the Cholesky decomposition expression associated to the underlying Delassus
    /// matrix.
    DelassusCholeskyExpression getDelassusCholeskyExpression() const
    {
      return DelassusCholeskyExpression(*this);
    }

    ///
    /// \brief Returns the Operational Space Inertia Matrix resulting from the decomposition.
    ///
    Matrix getOperationalSpaceInertiaMatrix() const
    {
      Matrix res(constraintDim(), constraintDim());
      getOperationalSpaceInertiaMatrix(res);
      return res;
    }

    template<typename MatrixType>
    void getOperationalSpaceInertiaMatrix(const Eigen::MatrixBase<MatrixType> & res_) const
    {
      MatrixType & res = PINOCCHIO_EIGEN_CONST_CAST(MatrixType, res_);
      //        typedef typename RowMatrix::ConstBlockXpr ConstBlockXpr;
      const auto U1 = U.topLeftCorner(constraintDim(), constraintDim())
                        .template triangularView<Eigen::UnitUpper>();

      const auto dim = constraintDim();
      typedef Eigen::Map<RowMatrix> MapRowMatrix;
      MapRowMatrix OSIMinv = MapRowMatrix(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, dim, dim));

      typedef Eigen::Map<Matrix> MapMatrix;
      MapMatrix U1inv = MapMatrix(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, dim, dim));

      PINOCCHIO_EIGEN_MALLOC_NOT_ALLOWED();
      U1inv.setIdentity();
      U1.solveInPlace(U1inv); // TODO: implement Sparse Inverse
      OSIMinv.noalias() = -U1inv.adjoint() * Dinv.head(dim).asDiagonal();
      res.noalias() = OSIMinv * U1inv;
      PINOCCHIO_EIGEN_MALLOC_ALLOWED();
    }

    Matrix getInverseMassMatrix() const
    {
      Matrix res(nv, nv);
      getInverseMassMatrix(res);
      return res;
    }

    template<typename MatrixType>
    void getInverseMassMatrix(const Eigen::MatrixBase<MatrixType> & res_) const
    {
      MatrixType & res = PINOCCHIO_EIGEN_CONST_CAST(MatrixType, res_);
      //        typedef typename RowMatrix::ConstBlockXpr ConstBlockXpr;
      const auto U4 = U.bottomRightCorner(nv, nv).template triangularView<Eigen::UnitUpper>();

      typedef Eigen::Map<RowMatrix> MapRowMatrix;
      MapRowMatrix Minv = MapRowMatrix(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, nv, nv));

      typedef Eigen::Map<Matrix> MapMatrix;
      MapMatrix U4inv = MapMatrix(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, nv, nv));

      PINOCCHIO_EIGEN_MALLOC_NOT_ALLOWED();
      U4inv.setIdentity();
      U4.solveInPlace(U4inv); // TODO: implement Sparse Inverse
      Minv.noalias() = U4inv.adjoint() * Dinv.tail(nv).asDiagonal();
      res.noalias() = Minv * U4inv;
      PINOCCHIO_EIGEN_MALLOC_ALLOWED();
    }

    template<typename MatrixType>
    void getJMinv(const Eigen::MatrixBase<MatrixType> & res_) const
    {
      PINOCCHIO_EIGEN_MALLOC_NOT_ALLOWED();
      MatrixType & res = PINOCCHIO_EIGEN_CONST_CAST(MatrixType, res_);
      const auto U4 = U.bottomRightCorner(nv, nv).template triangularView<Eigen::UnitUpper>();
      auto U2 = U.topRightCorner(constraintDim(), nv);

      typedef Eigen::Map<Matrix> MapMatrix;
      MapMatrix U4inv = MapMatrix(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, nv, nv));

      U4inv.setIdentity();
      U4.solveInPlace(U4inv); // TODO: implement Sparse Inverse
      res.noalias() = U2 * U4inv;
      PINOCCHIO_EIGEN_MALLOC_ALLOWED();
    }

    ///
    /// \brief Computes the Cholesky decompostion of the augmented matrix containing the KKT matrix
    ///        related to the system mass matrix and the Jacobians of the contact patches contained
    ///        in the vector of RigidConstraintModel named contact_models.
    ///
    /// \param[in] model Model of the dynamical system
    /// \param[in] data Data related to model containing the computed mass matrix and the Jacobian
    /// of the kinematic tree
    /// \param[in] contact_models Vector containing the contact models (which
    /// frame is in contact and the type of contact: ponctual, 6D rigid, etc.)
    /// \param[in,out] contact_datas Vector containing the contact data related to the
    /// contact_models.
    /// \param[in] mu Positive regularization factor allowing to enforce the definite property of
    /// the KKT matrix.
    ///
    /// \remarks The mass matrix and the Jacobians of the dynamical system should have been computed
    /// first. This can be achieved by simply calling pinocchio::crba.
    /// The `resize` method should have been called before calling this method if the size of the
    /// constraints changed.
    ///
    // TODO Remove when API is stabilized
    PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH
    PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS
    template<
      typename S1,
      int O1,
      template<typename, int> class JointCollectionTpl,
      class ConstraintModel,
      class ConstraintModelAllocator,
      class ConstraintData,
      class ConstraintDataAllocator>
    void compute(
      const ModelTpl<S1, O1, JointCollectionTpl> & model,
      DataTpl<S1, O1, JointCollectionTpl> & data,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & contact_models,
      std::vector<ConstraintData, ConstraintDataAllocator> & contact_datas,
      const S1 mu = S1(0.))
    {
      compute(model, data, contact_models, contact_datas, Vector::Constant(constraintDim(), mu));
    }

    ///
    /// \brief Computes the Cholesky decompostion of the augmented matrix containing the KKT matrix
    ///        related to the system mass matrix and the Jacobians of the contact patches contained
    ///        in the vector of ConstraintModel named contact_models.
    ///
    /// \param[in] model Model of the dynamical system
    /// \param[in] data Data related to model containing the computed mass matrix and the Jacobian
    /// of the kinematic tree
    /// \param[in] contact_models Vector containing the contact models (which
    /// frame is in contact and the type of contact: ponctual, 6D rigid, etc.)
    /// \param[in,out] contact_datas Vector containing the contact data related to the
    /// contact_models.
    /// \param[in] mu Positive regularization factor allowing to enforce the definite property of
    /// the KKT matrix.
    ///
    /// \remarks The mass matrix and the Jacobians of the dynamical system should have been computed
    /// first. This can be achieved by simply calling pinocchio::crba.
    /// The `resize` method should have been called before calling this method if the size of the
    /// constraints changed.
    ///
    template<
      typename S1,
      int O1,
      template<typename, int> class JointCollectionTpl,
      template<typename T> class Holder,
      class ConstraintModelAllocator,
      class ConstraintModel,
      class ConstraintDataAllocator,
      class ConstraintData>
    void compute(
      const ModelTpl<S1, O1, JointCollectionTpl> & model,
      DataTpl<S1, O1, JointCollectionTpl> & data,
      const std::vector<Holder<const ConstraintModel>, ConstraintModelAllocator> & contact_models,
      std::vector<Holder<ConstraintData>, ConstraintDataAllocator> & contact_datas,
      const S1 mu = S1(0.))
    {
      compute(model, data, contact_models, contact_datas, Vector::Constant(constraintDim(), mu));
    }
    PINOCCHIO_COMPILER_DIAGNOSTIC_POP

    ///
    /// \brief Computes the Cholesky decompostion of the augmented matrix containing the KKT matrix
    ///        related to the system mass matrix and the Jacobians of the contact patches contained
    ///        in the vector of ConstraintModel named contact_models.
    ///
    /// \param[in] model Model of the dynamical system
    /// \param[in] data Data related to model containing the computed mass matrix and the Jacobian
    /// of the kinematic tree
    /// \param[in] contact_models Vector containing the contact models (which
    /// frame is in contact and the type of contact: ponctual, 6D rigid, etc.)
    /// \param[in,out] contact_datas Vector containing the contact data related to the
    /// contact_models.
    /// \param[in] mu Positive regularization factor allowing to enforce the definite property of
    /// the KKT matrix.
    ///
    /// \remarks The mass matrix and the Jacobians of the dynamical system should have been computed
    /// first. This can be achieved by simply calling pinocchio::crba.
    /// The `resize` method should have been called before calling this method if the size of the
    /// constraints changed.
    ///
    template<
      typename S1,
      int O1,
      template<typename, int> class JointCollectionTpl,
      class ConstraintModel,
      class ConstraintModelAllocator,
      class ConstraintData,
      class ConstraintDataAllocator,
      typename VectorLike>
    void compute(
      const ModelTpl<S1, O1, JointCollectionTpl> & model,
      DataTpl<S1, O1, JointCollectionTpl> & data,
      const std::vector<ConstraintModel, ConstraintModelAllocator> & contact_models,
      std::vector<ConstraintData, ConstraintDataAllocator> & contact_datas,
      const Eigen::MatrixBase<VectorLike> & mus)
    {
      typedef std::reference_wrapper<const ConstraintModel> WrappedConstraintModelType;
      typedef std::vector<WrappedConstraintModelType> WrappedConstraintModelVector;

      WrappedConstraintModelVector wrapped_constraint_models(
        contact_models.cbegin(), contact_models.cend());

      typedef std::reference_wrapper<ConstraintData> WrappedConstraintDataType;
      typedef std::vector<WrappedConstraintDataType> WrappedConstraintDataVector;

      WrappedConstraintDataVector wrapped_constraint_datas(
        contact_datas.begin(), contact_datas.end());

      compute(model, data, wrapped_constraint_models, wrapped_constraint_datas, mus);
    }

    ///
    /// \brief Computes the Cholesky decompostion of the augmented matrix containing the KKT matrix
    ///        related to the system mass matrix and the Jacobians of the contact patches contained
    ///        in the vector of onstraintModel named contact_models.
    ///
    /// \param[in] model Model of the dynamical system
    /// \param[in] data Data related to model containing the computed mass matrix and the Jacobian
    /// of the kinematic tree
    /// \param[in] contact_models Vector containing the contact models (which
    /// frame is in contact and the type of contact: ponctual, 6D rigid, etc.)
    /// \param[in,out] contact_datas Vector containing the contact data related to the
    /// contact_models.
    /// \param[in] mu Positive regularization factor allowing to enforce the definite property of
    /// the KKT matrix.
    ///
    /// \remarks The mass matrix and the Jacobians of the dynamical system should have been computed
    /// first. This can be achieved by simply calling pinocchio::crba.
    /// The `resize` method should have been called before calling this method if the size of the
    /// constraints changed.
    ///
    template<
      typename S1,
      int O1,
      template<typename, int> class JointCollectionTpl,
      template<typename T> class Holder,
      class ConstraintModel,
      class ConstraintModelAllocator,
      class ConstraintData,
      class ConstraintDataAllocator,
      typename VectorLike>
    void compute(
      const ModelTpl<S1, O1, JointCollectionTpl> & model,
      DataTpl<S1, O1, JointCollectionTpl> & data,
      const std::vector<Holder<const ConstraintModel>, ConstraintModelAllocator> & contact_models,
      std::vector<Holder<ConstraintData>, ConstraintDataAllocator> & contact_datas,
      const Eigen::MatrixBase<VectorLike> & mus);

    ///
    /// \brief Update the compliance terms on the upper left block part of the KKT matrix. The
    /// compliance terms should be all positives.
    ///
    /// \param[in] compliance Vector of physical compliance for the constraints.
    ///
    template<typename VectorLike>
    void updateCompliance(const Eigen::MatrixBase<VectorLike> & compliance);

    ///
    /// \brief Update the compliance term on the upper left block part of the KKT matrix. The
    /// compliance terms should be all positives.
    ///
    /// \param[in] compliance The physical compliance for the constraints.
    ///
    void updateCompliance(const Scalar & compliance);

    ///
    /// \brief Returns the current compliance vector.
    ///
    const typename EigenStorageVector::MapType getCompliance() const
    {
      return compliance;
    }

    ///
    /// \brief Update the damping terms on the upper left block part of the KKT matrix. The damping
    /// terms should be all positives.
    ///
    /// \param[in] mus Vector of positive regularization factor allowing to enforce the definite
    /// property of the KKT matrix.
    ///
    template<typename VectorLike>
    void updateDamping(const Eigen::MatrixBase<VectorLike> & mus);

    ///
    /// \brief Update the damping term on the upper left block part of the KKT matrix. The damping
    /// terms should be all positives.
    ///
    /// \param[in] mu Regularization factor allowing to enforce the definite property of the KKT
    /// matrix.
    ///
    void updateDamping(const Scalar & mu);

    ///
    /// \brief Returns the current damping vector.
    ///
    const typename EigenStorageVector::MapType getDamping() const
    {
      return damping;
    }

    /// \brief Size of the decomposition
    Eigen::DenseIndex size() const
    {
      return D.size();
    }

    /// \brief Returns the total dimension of the constraints contained in the Cholesky
    /// factorization
    Eigen::DenseIndex constraintDim() const
    {
      return size() - nv;
    }

    ///
    ///  \brief Computes the solution of \f$ A x = b \f$ where *this is the Cholesky decomposition
    /// of A.         "in-place" version of ContactCholeskyDecompositionTpl::solve(b) where the
    /// result is written in b.
    ///        This functions takes as input the vector b, and returns the solution \f$ x = A^-1 b
    ///        \f$.
    ///
    /// \param[inout] mat The right-and-side term which also contains the solution of the linear
    /// system.
    ///
    /// \sa ContactCholeskyDecompositionTpl::solve
    template<typename MatrixLike>
    void solveInPlace(const Eigen::MatrixBase<MatrixLike> & mat) const;

    ///
    ///  \brief Computes the solution of \f$ A x = b \f$ where *this is the Cholesky decomposition
    /// of A.
    ///        This functions takes as input the vector b, and returns the solution \f$ x = A^-1 b
    ///        \f$.
    ///
    /// \param[inout] mat The right-and-side term.
    ///
    /// \sa ContactCholeskyDecompositionTpl::solveInPlace
    template<typename MatrixLike>
    Matrix solve(const Eigen::MatrixBase<MatrixLike> & mat) const;

    ///
    ///  \brief Retrieves the Cholesky decomposition of the Mass Matrix contained in *this.
    ///
    /// \param[in] model Model of the dynamical system.
    ///
    template<typename S1, int O1, template<typename, int> class JointCollectionTpl>
    ContactCholeskyDecompositionTpl
    getMassMatrixChoeslkyDecomposition(const ModelTpl<S1, O1, JointCollectionTpl> & model) const;

    ///@{
    /// \brief Vectorwize operations
    template<typename MatrixLike>
    void Uv(const Eigen::MatrixBase<MatrixLike> & mat) const;

    template<typename MatrixLike>
    void Utv(const Eigen::MatrixBase<MatrixLike> & mat) const;

    template<typename MatrixLike>
    void Uiv(const Eigen::MatrixBase<MatrixLike> & mat) const;

    template<typename MatrixLike>
    void Utiv(const Eigen::MatrixBase<MatrixLike> & mat) const;
    ///@}

    /// \brief Returns the matrix resulting from the decomposition
    Matrix matrix() const;

    /// \brief Fill the input matrix with the matrix resulting from the decomposition
    template<typename MatrixType>
    void matrix(const Eigen::MatrixBase<MatrixType> & res) const;

    /// \brief Returns the inverse matrix resulting from the decomposition
    Matrix inverse() const;

    /// \brief Fill the input matrix with the inverse matrix resulting from the decomposition
    template<typename MatrixType>
    void inverse(const Eigen::MatrixBase<MatrixType> & res) const;

    // data
    EigenStorageVector D_storage;
    typename EigenStorageVector::RefMapType D;
    EigenStorageVector Dinv_storage;
    typename EigenStorageVector::RefMapType Dinv;
    EigenStorageRowMatrix U_storage;
    typename EigenStorageRowMatrix::RefMapType U;

    ///@{
    /// \brief Friend algorithms
    template<typename MatrixLike, int ColsAtCompileTime>
    friend struct details::UvAlgo;

    template<typename MatrixLike, int ColsAtCompileTime>
    friend struct details::UtvAlgo;

    template<typename MatrixLike, int ColsAtCompileTime>
    friend struct details::UivAlgo;

    template<typename MatrixLike, int ColsAtCompileTime>
    friend struct details::UtivAlgo;

    // TODO Remove when API is stabilized
    PINOCCHIO_COMPILER_DIAGNOSTIC_PUSH
    PINOCCHIO_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS
    template<typename Scalar, int Options, typename VectorLike>
    friend VectorLike & details::inverseAlgo(
      const ContactCholeskyDecompositionTpl<Scalar, Options> & chol,
      const Eigen::DenseIndex col,
      const Eigen::MatrixBase<VectorLike> & vec);
    ///@}

    template<typename S1, int O1>
    bool operator==(const ContactCholeskyDecompositionTpl<S1, O1> & other) const;

    template<typename S1, int O1>
    bool operator!=(const ContactCholeskyDecompositionTpl<S1, O1> & other) const;
    PINOCCHIO_COMPILER_DIAGNOSTIC_POP

  protected:
    EigenIndexVector parents_fromRow;
    EigenIndexVector nv_subtree_fromRow;

    EigenStorageVector DUt_storage;
    typename EigenStorageVector::RefMapType DUt; // temporary containing the results of D * U^t

    /// \brief Dimension of the tangent of the configuration space of the model
    Eigen::DenseIndex nv;

    VectorOfSliceVector rowise_sparsity_pattern;

    /// \brief Store the current value of the physical compliance
    EigenStorageVector compliance_storage;
    typename EigenStorageVector::RefMapType compliance;

    /// \brief Store the current damping value
    EigenStorageVector damping_storage;
    typename EigenStorageVector::RefMapType damping;
  };

} // namespace pinocchio

// Because of a GCC bug we should NEVER define a function that use ContactCholeskyDecompositionTpl
// before doing the explicit template instantiation.
// If we don't take care, GCC will not accept any visibility attribute when declaring the
// explicit template instantiation of the ContactCholeskyDecompositionTpl class.
// The warning message will look like this: type attributes ignored after type is already defined
// [-Wattributes] A minimal code example is added on the PR
// (https://github.com/stack-of-tasks/pinocchio/pull/2469)
#if PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION
  #include "pinocchio/algorithm/contact-cholesky.txx"
#endif // PINOCCHIO_ENABLE_TEMPLATE_INSTANTIATION

#include "pinocchio/algorithm/contact-cholesky.hxx"
#include "pinocchio/algorithm/delassus-operator-cholesky-expression.hpp"

#endif // ifndef __pinocchio_algorithm_contact_cholesky_hpp__
