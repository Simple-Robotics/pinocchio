//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_delassus_operator_compliant_hpp__
#define __pinocchio_algorithm_delassus_operator_compliant_hpp__

#include "pinocchio/algorithm/fwd.hpp"
#include "pinocchio/math/arithmetic-operators.hpp"
#include "pinocchio/algorithm/delassus-operator-base.hpp"

namespace pinocchio
{

  template<typename DelassusOperator, typename ComplianceType>
  struct DelassusOperatorCompliantTpl;

  template<typename DelassusOperator, typename ComplianceType>
  struct traits<DelassusOperatorCompliantTpl<DelassusOperator, ComplianceType>>
  : traits<DelassusOperator>
  {
  };

  template<typename DelassusOperator, typename ComplianceType>
  struct DelassusOperatorCompliantTpl
  : DelassusOperatorBase<DelassusOperatorCompliantTpl<DelassusOperator, ComplianceType>>
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef DelassusOperatorCompliantTpl Self;
    typedef DelassusOperatorBase<Self> Base;

    typedef typename traits<Self>::Matrix Matrix;
    typedef typename traits<Self>::Vector Vector;
    typedef typename traits<Self>::Scalar Scalar;

    DelassusOperatorCompliantTpl(
      DelassusOperatorBase<DelassusOperator> & delassus, ComplianceType & compliance)
    : m_delassus(delassus.derived())
    , m_compliance(compliance)
    , m_tmp_vec(compliance.cols())
    {
      PINOCCHIO_CHECK_ARGUMENT_SIZE(m_delassus.cols(), m_compliance.rows());
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
    void updateCompliance(const Eigen::MatrixBase<VectorLike> & vec)
    {
      m_compliance = vec;
    }

    void updateCompliance(const Scalar mu)
    {
      this->updateCompliance(Vector::Constant(ref().size(), mu));
    }

    template<typename VectorLike>
    void updateDamping(const Eigen::MatrixBase<VectorLike> & vec)
    {
      // G + R + mu * Id
      m_tmp_vec = m_compliance + vec;
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
      ref().solveInPlace(mat_);
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
      ref().applyOnTheRight(res_, x);
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

    Matrix matrix(bool enforce_symmetry = false) const
    {
      return m_delassus.matrix(enforce_symmetry) + m_compliance.asDiagonal().toDenseMatrix();
    }

    const Vector & getDamping() const
    {
      m_tmp_vec = ref().getDamping() - m_compliance;
      return m_tmp_vec;
    }

    const Vector & getCompliance() const
    {
      return m_compliance;
    }

  protected:
    DelassusOperator & m_delassus;
    ComplianceType & m_compliance;
    Vector m_tmp_vec;

  }; // struct DelassusOperatorCompliant

} // namespace pinocchio

#endif // ifndef __pinocchio_algorithm_delassus_operator_compliant_hpp__
