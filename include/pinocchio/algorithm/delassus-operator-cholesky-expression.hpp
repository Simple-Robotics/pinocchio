//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_delassus_operator_cholesky_expression_hpp__
#define __pinocchio_algorithm_delassus_operator_cholesky_expression_hpp__

#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/algorithm/delassus-operator-base.hpp"
#include "pinocchio/algorithm/contact-cholesky.hpp"

namespace pinocchio
{
  template<typename ContactCholeskyDecomposition>
  struct traits<DelassusCholeskyExpressionTpl<ContactCholeskyDecomposition>>
  {
    enum
    {
      RowsAtCompileTime = Eigen::Dynamic
    };
    typedef typename ContactCholeskyDecomposition::Scalar Scalar;
    typedef typename ContactCholeskyDecomposition::Matrix Matrix;
    typedef typename ContactCholeskyDecomposition::Vector Vector;
  };

  template<typename _ContactCholeskyDecomposition>
  struct DelassusCholeskyExpressionTpl
  : DelassusOperatorBase<DelassusCholeskyExpressionTpl<_ContactCholeskyDecomposition>>
  {
    typedef _ContactCholeskyDecomposition ContactCholeskyDecomposition;
    typedef typename ContactCholeskyDecomposition::Scalar Scalar;
    typedef typename ContactCholeskyDecomposition::Vector Vector;
    typedef typename ContactCholeskyDecomposition::Matrix Matrix;
    typedef typename ContactCholeskyDecomposition::RowMatrix RowMatrix;
    typedef DelassusCholeskyExpressionTpl<_ContactCholeskyDecomposition> Self;
    typedef DelassusOperatorBase<Self> Base;

    typedef
      typename SizeDepType<Eigen::Dynamic>::template BlockReturn<RowMatrix>::Type RowMatrixBlockXpr;
    typedef typename SizeDepType<Eigen::Dynamic>::template BlockReturn<RowMatrix>::ConstType
      RowMatrixConstBlockXpr;

    enum
    {
      RowsAtCompileTime = traits<DelassusCholeskyExpressionTpl>::RowsAtCompileTime
    };

    explicit DelassusCholeskyExpressionTpl(const ContactCholeskyDecomposition & self)
    : Base()
    , self(self)
    {
    }

    template<typename MatrixIn, typename MatrixOut>
    void applyOnTheRight(
      const Eigen::MatrixBase<MatrixIn> & x, const Eigen::MatrixBase<MatrixOut> & res) const
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(x.rows(), self.constraintDim());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(res.rows(), self.constraintDim());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(res.cols(), x.cols());

      const auto U1 = self.U.topLeftCorner(self.constraintDim(), self.constraintDim());
      {
        PINOCCHIO_EIGEN_MALLOC_NOT_ALLOWED();
        typedef Eigen::Map<RowMatrix> MapType;
        MapType tmp_mat = MapType(PINOCCHIO_EIGEN_MAP_ALLOCA(Scalar, x.rows(), x.cols()));
        //            tmp_mat.noalias() = U1.adjoint() * x;
        triangularMatrixMatrixProduct<Eigen::UnitLower>(U1.adjoint(), x.derived(), tmp_mat);

        // The following commented lines produced some memory allocation.
        // Should be replaced by a manual loop
        //          tmp_mat.array().colwise() *= -self.D.head(self.constraintDim()).array();
        for (Eigen::DenseIndex i = 0; i < x.cols(); ++i)
          tmp_mat.col(i).array() *= -self.D.head(self.constraintDim()).array();

        //            res.const_cast_derived().noalias() = U1 * tmp_mat;
        triangularMatrixMatrixProduct<Eigen::UnitUpper>(U1, tmp_mat, res.const_cast_derived());
        PINOCCHIO_EIGEN_MALLOC_ALLOWED();
      }
    }

    template<typename MatrixDerived>
    void solveInPlace(const Eigen::MatrixBase<MatrixDerived> & x) const
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(x.rows(), self.constraintDim());

      const auto U1 = self.U.topLeftCorner(self.constraintDim(), self.constraintDim())
                        .template triangularView<Eigen::UnitUpper>();

      PINOCCHIO_EIGEN_MALLOC_NOT_ALLOWED();
      U1.solveInPlace(x.const_cast_derived());

      // The following commented lines produced some memory allocation.
      // Should be replaced by a manual loop
      //        x.const_cast_derived().array().colwise() *=
      //        -self.Dinv.head(self.constraintDim()).array();
      for (Eigen::DenseIndex i = 0; i < x.cols(); ++i)
        x.const_cast_derived().col(i).array() *= -self.Dinv.head(self.constraintDim()).array();

      U1.adjoint().solveInPlace(x);
      PINOCCHIO_EIGEN_MALLOC_ALLOWED();
    }

    template<typename MatrixDerivedIn, typename MatrixDerivedOut>
    void solve(
      const Eigen::MatrixBase<MatrixDerivedIn> & x,
      const Eigen::MatrixBase<MatrixDerivedOut> & res) const
    {
      res.const_cast_derived() = x;
      solveInPlace(res.const_cast_derived());
    }

    template<typename MatrixDerived>
    typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixDerived)
      solve(const Eigen::MatrixBase<MatrixDerived> & x) const
    {
      typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixDerived) ReturnType;
      ReturnType res(self.constraintDim(), x.cols());
      solve(x.derived(), res);
      return res;
    }

    /// \brief Returns the Constraint Cholesky decomposition associated to this
    /// DelassusCholeskyExpression.
    const ContactCholeskyDecomposition & cholesky() const
    {
      return self;
    }

    Matrix matrix() const
    {
      return self.getInverseOperationalSpaceInertiaMatrix();
    }

    /// \brief Fill the input matrix with the matrix resulting from the decomposition
    template<typename MatrixType>
    void matrix(const Eigen::MatrixBase<MatrixType> & mat) const
    {
      return self.getInverseOperationalSpaceInertiaMatrix(mat.const_cast_derived());
    }

    ///
    /// \brief Returns the current damping vector.
    ///
    const Vector & getDamping() const
    {
      return self.getDamping();
    }

    Matrix inverse() const
    {
      return self.getOperationalSpaceInertiaMatrix();
    }

    ///
    /// \brief Add a damping term to the diagonal of the Delassus matrix. The damping terms should
    /// be all positives.
    ///
    /// \param[in] mus Vector of positive regularization factor allowing to enforce the definite
    /// positiveness of the matrix.
    ///
    template<typename VectorLike>
    void updateDamping(const Eigen::MatrixBase<VectorLike> & mus)
    {
      const_cast<ContactCholeskyDecomposition &>(self).updateDamping(mus);
    }

    ///
    /// \brief Add a damping term to the diagonal of the Delassus matrix. The damping term should be
    /// positive.
    ///
    /// \param[in] mu Regularization factor allowing to enforce the definite positiveness of the
    /// matrix.
    ///
    void updateDamping(const Scalar & mu)
    {
      const_cast<ContactCholeskyDecomposition &>(self).updateDamping(mu);
    }

    Eigen::DenseIndex size() const
    {
      return self.constraintDim();
    }
    Eigen::DenseIndex rows() const
    {
      return size();
    }
    Eigen::DenseIndex cols() const
    {
      return size();
    }

  protected:
    const ContactCholeskyDecomposition & self;
  }; // DelassusCholeskyExpression

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_delassus_operator_cholesky_expression_hpp__
