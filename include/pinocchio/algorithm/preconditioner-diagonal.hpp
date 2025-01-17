//
// Copyright (c) 2025 INRIA
//

#ifndef __pinocchio_algorithm_preconditioner_diagonal_hpp__
#define __pinocchio_algorithm_preconditioner_diagonal_hpp__

#include "pinocchio/algorithm/preconditioner-base.hpp"
namespace pinocchio
{

  template<typename PreconditionerVectorLike>
  struct PreconditionerDiagonal
  : PreconditionerBase<PreconditionerDiagonal<PreconditionerVectorLike>>
  {

    explicit PreconditionerDiagonal(const PreconditionerVectorLike & diagonal)
    : m_preconditioner_diagonal(diagonal)
    , m_preconditioner_square(diagonal)
    {
      m_preconditioner_square.array() *= diagonal.array();
    }

    template<typename MatrixIn, typename MatrixOut>
    void
    scale(const Eigen::MatrixBase<MatrixIn> & x, const Eigen::MatrixBase<MatrixOut> & res) const
    {
      auto & res_ = res.const_cast_derived();
      res_.array() = x.array() / m_preconditioner_diagonal.array();
    }

    template<typename MatrixIn>
    void scaleInPlace(const Eigen::MatrixBase<MatrixIn> & x) const
    {
      auto & x_ = x.const_cast_derived();
      x_.array() = x.array() / m_preconditioner_diagonal.array();
    }

    template<typename MatrixIn, typename MatrixOut>
    void scaleSquare(
      const Eigen::MatrixBase<MatrixIn> & x, const Eigen::MatrixBase<MatrixOut> & res) const
    {
      auto & res_ = res.const_cast_derived();
      res_.array() = x.array() / m_preconditioner_square.array();
    }

    template<typename MatrixIn, typename MatrixOut>
    void
    unscale(const Eigen::MatrixBase<MatrixIn> & x, const Eigen::MatrixBase<MatrixOut> & res) const
    {
      auto & res_ = res.const_cast_derived();
      res_.array() = x.array() * m_preconditioner_diagonal.array();
    }

    template<typename MatrixIn>
    void unscaleInPlace(const Eigen::MatrixBase<MatrixIn> & x) const
    {
      auto & x_ = x.const_cast_derived();
      x_.array() *= m_preconditioner_diagonal.array();
    }

    template<typename MatrixIn, typename MatrixOut>
    void unscaleSquare(
      const Eigen::MatrixBase<MatrixIn> & x, const Eigen::MatrixBase<MatrixOut> & res) const
    {
      auto & res_ = res.const_cast_derived();
      res_.array() = x.array() * m_preconditioner_square.array();
    }

    Eigen::DenseIndex rows() const
    {
      return m_preconditioner_diagonal.size();
    }
    Eigen::DenseIndex cols() const
    {
      return m_preconditioner_diagonal.size();
    }

    void setDiagonal(const PreconditionerVectorLike & x)
    {
      m_preconditioner_diagonal = x;
      m_preconditioner_square.array() = x.array() * x.array();
    }

    const PreconditionerVectorLike & getDiagonal() const
    {
      return m_preconditioner_diagonal;
    }

  protected:
    PreconditionerVectorLike m_preconditioner_diagonal;
    PreconditionerVectorLike m_preconditioner_square;

  }; // struct PreconditionerDiagonal

} // namespace pinocchio

#endif // #ifndef __pinocchio_algorithm_preconditioner_diagonal_hpp__
