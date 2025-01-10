//
// Copyright (c) 2024 INRIA
//

#ifndef __pinocchio_algorithm_delassus_operator_preconditionned_hpp__
#define __pinocchio_algorithm_delassus_operator_preconditionned_hpp__

#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/math/eigenvalues.hpp"
#include "pinocchio/math/arithmetic-operators.hpp"
#include "pinocchio/algorithm/delassus-operator-base.hpp"

namespace pinocchio
{

  template<typename DelassusOperatorDerived>
  struct DelassusOperatorPreconditionnedTpl : DelassusOperatorBase<DelassusOperatorDerived>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef DelassusOperatorPreconditionnedTpl Self;
    typedef DelassusOperatorBase<DelassusOperatorDerived> Base;
    typedef typename traits<DelassusOperatorDerived>::Vector Vector;
    typedef typename traits<DelassusOperatorDerived>::Scalar Scalar;

    template<typename VectorDerived>
    explicit DelassusOperatorPreconditionnedTpl(
      Base & delassus, const Eigen::MatrixBase<VectorDerived> & preconditioner)
    : m_preconditionner(preconditioner)
    , m_damping(preconditioner)
    , m_delassus(delassus)
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(m_delassus.rows(), m_preconditionner.size());
    }

    DelassusOperatorDerived & ref()
    {
      return m_delassus.derived();
    }
    const DelassusOperatorDerived & ref() const
    {
      return m_delassus.derived();
    }

    template<typename VectorLike>
    void updateDamping(const Eigen::MatrixBase<VectorLike> & vec)
    {
      // G_bar + mu * Id = P * (G + mu * P^{-2}) * P
      m_damping.array() = vec.array() / m_preconditionner.array();
      m_damping.array() /= m_preconditionner.array();
      ref().updateDamping(m_damping);
    }

    void updateDamping(const Scalar mu)
    {
      this->updateDamping(Vector::Constant(ref().size(), mu));
    }

    template<typename MatrixLike>
    void solveInPlace(const Eigen::MatrixBase<MatrixLike> & mat) const
    {
      auto mat_ = mat.const_cast_derived();
      mat_.array() /= m_preconditionner.array();
      ref().solveInPlace(mat_.const_cast_derived());
      mat_.array() /= m_preconditionner.array();
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
      auto res_ = res.const_cast_derived();
      res_.array() *= m_preconditionner.array();
      ref().applyOnTheRight(x, res_);
      res_.array() *= m_preconditionner.array();
    }

  protected:
    Base & m_delassus;
    const Vector & m_preconditionner;
    Vector m_damping;

  }; // struct DelassusOperatorPreconditionned

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_delassus_operator_preconditionned_hpp__
