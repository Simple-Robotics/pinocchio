//
// Copyright (c) 2024-2025 INRIA
//

#ifndef __pinocchio_algorithm_delassus_operator_dense_hpp__
#define __pinocchio_algorithm_delassus_operator_dense_hpp__

#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/algorithm/delassus-operator-base.hpp"
#include "pinocchio/algorithm/contact-cholesky.hpp"

namespace pinocchio
{

  template<typename _Scalar, int _Options>
  struct traits<DelassusOperatorDenseTpl<_Scalar, _Options>>
  {
    typedef _Scalar Scalar;
    enum
    {
      Options = _Options,
      RowsAtCompileTime = Eigen::Dynamic
    };

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> Matrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Options> Vector;

    typedef const Vector & getDampingReturnType;
  };

  template<typename _Scalar, int _Options>
  struct DelassusOperatorDenseTpl
  : DelassusOperatorBase<DelassusOperatorDenseTpl<_Scalar, _Options>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef _Scalar Scalar;
    typedef DelassusOperatorDenseTpl Self;
    enum
    {
      Options = _Options,
      RowsAtCompileTime = traits<DelassusOperatorDenseTpl>::RowsAtCompileTime
    };

    typedef typename traits<Self>::Matrix Matrix;
    typedef typename traits<Self>::Vector Vector;
    typedef Eigen::LLT<Matrix> CholeskyDecomposition;
    typedef DelassusOperatorBase<Self> Base;

    template<typename MatrixDerived>
    explicit DelassusOperatorDenseTpl(const Eigen::MatrixBase<MatrixDerived> & mat)
    : Base()
    , delassus_matrix(mat)
    , mat_tmp(mat.rows(), mat.cols())
    , llt(mat)
    , compliance(Vector::Zero(mat.rows()))
    , damping(Vector::Zero(mat.rows()))
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(mat.rows(), mat.cols());
    }

    template<typename ContactCholeskyDecomposition>
    explicit DelassusOperatorDenseTpl(
      const DelassusCholeskyExpressionTpl<ContactCholeskyDecomposition> & delassus_expression,
      const bool enforce_symmetry = false)
    : Base()
    , delassus_matrix(delassus_expression.matrix(enforce_symmetry))
    , mat_tmp(delassus_expression.rows(), delassus_expression.cols())
    , llt(delassus_matrix)
    , compliance(delassus_expression.getCompliance())
    , damping(delassus_expression.getDamping())
    {
      delassus_matrix -= delassus_expression.getDamping().asDiagonal();
      delassus_matrix -= delassus_expression.getCompliance().asDiagonal();
    }

    template<typename VectorLike>
    void updateCompliance(const Eigen::MatrixBase<VectorLike> & vec)
    {
      compliance = vec;
      updateDamping(getDamping());
    }

    void updateCompliance(const Scalar & compliance)
    {
      updateCompliance(Vector::Constant(size(), compliance));
    }

    template<typename VectorLike>
    void updateDamping(const Eigen::MatrixBase<VectorLike> & vec)
    {
      damping = vec;
      mat_tmp = delassus_matrix;
      mat_tmp += vec.asDiagonal();
      mat_tmp += compliance.asDiagonal();
      llt.compute(mat_tmp);
    }

    void updateDamping(const Scalar & mu)
    {
      updateDamping(Vector::Constant(size(), mu));
    }

    template<typename MatrixLike>
    void solveInPlace(const Eigen::MatrixBase<MatrixLike> & mat) const
    {
      llt.solveInPlace(mat.const_cast_derived());
    }

    template<typename MatrixLike>
    typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixLike)
      solve(const Eigen::MatrixBase<MatrixLike> & mat) const
    {
      typename PINOCCHIO_EIGEN_PLAIN_TYPE(MatrixLike) res(mat);
      solveInPlace(res);
      return res;
    }

    template<typename MatrixDerivedIn, typename MatrixDerivedOut>
    void solve(
      const Eigen::MatrixBase<MatrixDerivedIn> & x,
      const Eigen::MatrixBase<MatrixDerivedOut> & res) const
    {
      res.const_cast_derived() = x;
      solveInPlace(res.const_cast_derived());
    }

    template<typename MatrixIn, typename MatrixOut>
    void applyOnTheRight(
      const Eigen::MatrixBase<MatrixIn> & x, const Eigen::MatrixBase<MatrixOut> & res_) const
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(x.rows(), size());
      MatrixOut & res = res_.const_cast_derived();
      res.noalias() = delassus_matrix * x;
      res.array() += damping.array() * x.array();
      res.array() += compliance.array() * x.array();
    }

    Eigen::DenseIndex size() const
    {
      return delassus_matrix.rows();
    }
    Eigen::DenseIndex rows() const
    {
      return delassus_matrix.rows();
    }
    Eigen::DenseIndex cols() const
    {
      return delassus_matrix.cols();
    }

    Matrix matrix(bool enforce_symmetry = false) const
    {
      mat_tmp = delassus_matrix;
      mat_tmp += damping.asDiagonal();
      mat_tmp += compliance.asDiagonal();
      if (enforce_symmetry)
      {
        enforceSymmetry(mat_tmp);
      }
      return mat_tmp;
    }

    Matrix undampedMatrix(bool enforce_symmetry = false) const
    {
      mat_tmp = delassus_matrix;
      mat_tmp += compliance.asDiagonal();
      if (enforce_symmetry)
      {
        enforceSymmetry(mat_tmp);
      }
      return mat_tmp;
    }

    const Vector & getCompliance() const
    {
      return compliance;
    }

    const Vector & getDamping() const
    {
      return damping;
    }

    Matrix inverse() const
    {
      Matrix res = Matrix::Identity(size(), size());
      llt.solveInPlace(res);
      return res;
    }

    bool operator==(const Self & other) const
    {
      if (&other == this)
        return true;
      return delassus_matrix == other.delassus_matrix && damping == other.damping
             && compliance == other.compliance;
    }

    bool operator!=(const Self & other) const
    {
      return !(*this == other);
    }

  protected:
    Matrix delassus_matrix;
    mutable Matrix mat_tmp;
    CholeskyDecomposition llt;
    Vector damping;
    Vector compliance;

  }; // struct DelassusOperatorDenseTpl

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_delassus_operator_dense_hpp__
