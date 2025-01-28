//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_delassus_operator_preconditioned_hpp__
#define __pinocchio_algorithm_delassus_operator_preconditioned_hpp__

#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/math/arithmetic-operators.hpp"
#include "pinocchio/algorithm/delassus-operator-base.hpp"
#include "pinocchio/algorithm/preconditioner-base.hpp"

namespace pinocchio
{

  template<typename DelassusOperator, typename PreconditionerType>
  struct DelassusOperatorPreconditionedTpl;

  template<typename DelassusOperator, typename PreconditionerType>
  struct traits<DelassusOperatorPreconditionedTpl<DelassusOperator, PreconditionerType>>
  : traits<DelassusOperator>
  {
  };

  template<typename DelassusOperator, typename PreconditionerType>
  struct DelassusOperatorPreconditionedTpl
  : DelassusOperatorBase<DelassusOperatorPreconditionedTpl<DelassusOperator, PreconditionerType>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef DelassusOperatorPreconditionedTpl Self;
    typedef DelassusOperatorBase<Self> Base;

    typedef typename traits<Self>::Matrix Matrix;
    typedef typename traits<Self>::Vector Vector;
    typedef typename traits<Self>::Scalar Scalar;

    DelassusOperatorPreconditionedTpl(
      DelassusOperatorBase<DelassusOperator> & delassus, const PreconditionerType & preconditioner)
    : m_delassus(delassus.derived())
    , m_preconditioner(preconditioner)
    , m_tmp_vec(preconditioner.cols())
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(m_delassus.rows(), m_preconditioner.cols());
      PINOCCHIO_CHECK_ARGUMENT_SIZE(m_delassus.cols(), m_preconditioner.rows());
    }

    DelassusOperator & ref()
    {
      return m_delassus;
    }
    const DelassusOperator & ref() const
    {
      return m_delassus;
    }

    template<typename VectorLike>
    void updateDamping(const Eigen::MatrixBase<VectorLike> & vec)
    {
      // G_bar + mu * Id = P * (G + mu * P^{-2}) * P
      m_preconditioner.scaleSquare(vec, m_tmp_vec);
      ref().updateDamping(m_tmp_vec);
    }

    void updateDamping(const Scalar mu)
    {
      this->updateDamping(Vector::Constant(ref().size(), mu));
    }

    template<typename MatrixLike>
    void solveInPlace(const Eigen::MatrixBase<MatrixLike> & mat) const
    {
      auto & mat_ = mat.const_cast_derived();
      m_preconditioner.scaleInPlace(mat_);
      ref().solveInPlace(mat_);
      m_preconditioner.scaleInPlace(mat_);
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
      const Eigen::MatrixBase<MatrixIn> & x, const Eigen::MatrixBase<MatrixOut> & res) const
    {
      auto & res_ = res.const_cast_derived();
      m_preconditioner.unscale(x, res_);
      ref().applyOnTheRight(res_, m_tmp_vec);
      m_preconditioner.unscale(m_tmp_vec, res_);
    }

    Eigen::DenseIndex size() const
    {
      return ref().size();
    }
    Eigen::DenseIndex rows() const
    {
      return ref().rows();
    }
    Eigen::DenseIndex cols() const
    {
      return ref().cols();
    }

    Matrix matrix() const
    {
      return m_preconditioner.getDiagonal().asDiagonal() * m_delassus.matrix()
             * m_preconditioner.getDiagonal().asDiagonal();
    }

    const Vector & getDamping() const
    {
      m_preconditioner.unscaleSquare(ref().getDamping(), m_tmp_vec);
      return m_tmp_vec;
    }

  protected:
    DelassusOperator & m_delassus;
    const PreconditionerType & m_preconditioner;
    Vector m_tmp_vec;

  }; // struct DelassusOperatorPreconditioned

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_delassus_operator_preconditioned_hpp__
